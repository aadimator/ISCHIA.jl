"""
Calculate significant co-occurring Ligand-Receptor pairs.

This function calculates co-occurring Ligand-Receptor (LR) pairs that are statistically significant based on expression levels and correlations in a spatial dataset.

Parameters:
- `adata::AnnData`: The (spatial) anndata dataset containing expression data.
- `COI::Vector{String}`: Cluster of Interest, a subset of spots to focus on.
- `Condition::Vector{String}`: Condition of interest within the dataset.
- `LR_list::Vector{String}`: List of ligands and receptors to consider.
- `LR_pairs::Vector{String}`: List of LR pairs to analyze.
- `exp_th::Real`: Expression threshold for binarizing the expression matrix.
- `corr_th::Real`: Correlation threshold for LR pairs.

Returns:
A dictionary containing:
- `"enriched_LRs"`: DataFrame of enriched LR pairs.
- `"cooccurrence_table"`: Co-occurrence analysis results.

"""
function find_enriched_LR_pairs(
    adata::AnnData, COI::Vector{String}, Condition::Vector{String},
    LR_list::Vector{String}, LR_pairs::Vector{String},
    exp_th::Real, corr_th::Real; cc_column::String="CompositionCluster_CC", condition_column::String="orig.ident")

    println("Preparing L-R presence/absence matrix")

    # Subset the expression matrix for the interested ligands and receptors
    spatial_obj_exp_LR_subset_raw = adata[:, LR_list]

    # Binarize the expression matrix based on the expression threshold
    spatial_obj_exp_LR_subset_raw_binary = spatial_obj_exp_LR_subset_raw.layers["counts"] .> exp_th
    spatial_obj_exp_LR_subset_raw.layers["binary"] = spatial_obj_exp_LR_subset_raw_binary

    LR_subset_raw_binary_mask_col = vec(sum(spatial_obj_exp_LR_subset_raw_binary, dims=1) .> 0)
    LR_subset_raw_binary_mask_row = vec(sum(spatial_obj_exp_LR_subset_raw_binary, dims=2) .> 0)

    LR_presence_absence = spatial_obj_exp_LR_subset_raw[LR_subset_raw_binary_mask_row, LR_subset_raw_binary_mask_col]

    # Filter spots based on COI and Condition
    mask = (adata.obs[:, cc_column] .∈ Ref(COI)) .& (adata.obs[:, condition_column] .∈ Ref(Condition))
    COI_spots = adata.obs_names[mask]
    rest_of_spots = setdiff(adata.obs_names, COI_spots)

    println("Calculating L-R pairs correlation")
    COI_cors_adata = spatial_obj_exp_LR_subset_raw[mask, :]
    COI_cors = cor(Array(COI_cors_adata.layers["counts"]))
    COI_cors[isnan.(COI_cors)] .= 0.0

    println("Preparing for cooccurrence")
    common_spots = intersect(LR_presence_absence.obs_names, COI_spots)
    coocur_COI = LR_presence_absence[common_spots, :]
    coocur_COI_exp = DataFrame(Matrix(transpose(coocur_COI.layers["binary"])), common_spots)

    println("Cooccurrence calculation starts...")
    cooccur_COI_res = calculate_cooccurrence_stats(Matrix(coocur_COI_exp), coocur_COI.var.name; spp_names=true)
    println("Cooccurrence calculation ended")

    println("\nSummary of cooccurrence results:")
    summarize_cooccur(cooccur_COI_res)

    # println("\nProbability table of cooccurrence results:")
    flush(stdout)
    # display(cooccur_COI_res.results)

    cooccur_res_df = cooccur_COI_res.results
    # Add a 'pair' column to the result DataFrame
    cooccur_res_df[!, :pair] = string.(cooccur_res_df.sp1_name, "_", cooccur_res_df.sp2_name)
    common_pairs = intersect(LR_pairs, cooccur_res_df.pair)

    COI_enriched_LRs = DataFrame(from=String[], to=String[], correlation=Float64[], ligand_FC=Float64[], Receptor_FC=Float64[])
    @showprogress "Calculate Significantly occurring pairs" for pair in common_pairs

        # Split the LR pair into individual ligand and receptor
        LR_pair_words = split(pair, "_")
        LR_pair_ligand = String(LR_pair_words[1])
        LR_pair_Receptor = String(LR_pair_words[2])

        # Mean expression of the ligand in the Cluster of Interest (COI) spots and rest of the spots
        ligand_exp_COI_mean = mean(adata[COI_spots, LR_pair_ligand].X)
        ligand_exp_otherspots_mean = mean(adata[rest_of_spots, LR_pair_ligand].X)
        # Calculate the ligand fold change (FC) by dividing COI mean by rest of the spots mean
        ligand_FC = round(ligand_exp_COI_mean / ligand_exp_otherspots_mean, digits=4)

        Receptor_exp_COI_mean = mean(adata[COI_spots, LR_pair_Receptor].X)
        Receptor_exp_otherspots_mean = mean(adata[rest_of_spots, LR_pair_Receptor].X)
        Receptor_FC = round(Receptor_exp_COI_mean / Receptor_exp_otherspots_mean, digits=4)

        # Retrieve the p-value for the pair from the co-occurrence results DataFrame
        pair_p = cooccur_res_df[cooccur_res_df.pair .== pair, :p_gt][1]

        # Find the indices of the ligand and receptor in the COI correlation matrix
        ligand_index = findfirst(==(LR_pair_ligand), COI_cors_adata.var_names)
        receptor_index = findfirst(==(LR_pair_Receptor), COI_cors_adata.var_names)

        # Check if the pair is significant (p-value < 0.05) and the correlation is above the threshold
        if pair_p < 0.05 && COI_cors[ligand_index, receptor_index] > corr_th
            added_row = DataFrame(from=[LR_pair_ligand], to=[LR_pair_Receptor], correlation=[COI_cors[ligand_index, receptor_index]], ligand_FC=[ligand_FC], Receptor_FC=[Receptor_FC])
            append!(COI_enriched_LRs, added_row)
        end
    end

    # Sort the enriched LRs by correlation in decreasing order
    sort!(COI_enriched_LRs, rev=true, [:correlation])

    # Add a 'pair' column to the enriched LRs DataFrame
    COI_enriched_LRs[!, :pair] = string.(COI_enriched_LRs[!, :from], "_", COI_enriched_LRs[!, :to])

    Output_dict = Dict("enriched_LRs" => COI_enriched_LRs, "cooccurrence_table" => cooccur_COI_res)
    return Output_dict
end


"""
Find LR (Ligand Receptor) pairs that are significantly co-occurring in one group and not in the other group.

# Arguments
- `group1_results`: Results from the EnrichedLRs function for Group 1.
- `group2_results`: Results from the EnrichedLRs function for Group 2.
- `group1_max_pval`: Maximum p-value threshold for significance levels of co-occurring LR pairs in Group 1.
- `group2_min_pval`: Minimum p-value threshold for non-significance levels of co-occurring LR pairs in Group 2.

# Returns
List of LR pairs enriched in Group 1 and not in Group 2.

# Example
```julia
result = find_differentially_cooccurring_LR_pairs(results_group1, results_group2, 0.05, 0.1)
"""
function find_differentially_cooccurring_LR_pairs(cooc_df_1, cooc_df_2, group1_max_pval, group2_min_pval)

    enriched_LR_pairs_group1 = DataFrame(pair=String[], group1_pval=Real[], group2_pval=Real[], pval_difference=Real[], observed_cooc=Int[])

    @showprogress for row in eachrow(cooc_df_1)
        if row.pair in cooc_df_2.pair
            group2_row = filter(r -> r.pair == row.pair, cooc_df_2)
            group1_pval = row.p_gt
            group2_pval = group2_row.p_gt[1]
            group1_observed_cooc = row.obs_cooccur
            group2_observed_cooc = group2_row.obs_cooccur[1]
            group1_expected_cooc = row.exp_cooccur
            group2_expected_cooc = group2_row.exp_cooccur[1]
            pval_difference = group2_pval - group1_pval

            if group1_pval < group1_max_pval && group2_pval > group2_min_pval &&
               group1_observed_cooc > 10 &&
               group1_observed_cooc != group1_expected_cooc && group2_observed_cooc != group2_expected_cooc &&
               group2_observed_cooc < group1_observed_cooc
                pair_data = DataFrame(pair=row.pair, group1_pval=group1_pval, group2_pval=group2_pval,
                    pval_difference=pval_difference, observed_cooc=group1_observed_cooc
                )
                append!(enriched_LR_pairs_group1, pair_data)
            end
        end
    end

    enriched_LR_pairs_group1_sorted = sort(enriched_LR_pairs_group1, :observed_cooc, rev=true)

    return enriched_LR_pairs_group1_sorted
end