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
function enriched_LRs(
    adata::AnnData, COI::Vector{String}, Condition::Vector{String},
    LR_list::Vector{String}, LR_pairs::Vector{String},
    exp_th::Real, corr_th::Real)

    println("Preparing L-R presence/absence matrix")

    # Extract the expression matrix from spatial_object
    spatial_object_exp = adata.layers["counts"]
    spatial_object_exp_norm = adata.X

    # Subset the expression matrix for the interested ligands and receptors
    spatial_obj_exp_LR_subset_raw = adata[:, in.(adata.var.name, Ref(LR_list))]

    # Binarize the expression matrix based on the expression threshold
    spatial_obj_exp_LR_subset_raw_binary = spatial_obj_exp_LR_subset_raw.layers["counts"] .> exp_th
    spatial_obj_exp_LR_subset_raw.layers["binary"] = spatial_obj_exp_LR_subset_raw_binary

    LR_subset_raw_binary_mask_col = vec(sum(spatial_obj_exp_LR_subset_raw_binary, dims=1) .> 0)
    LR_subset_raw_binary_mask_row = vec(sum(spatial_obj_exp_LR_subset_raw_binary, dims=2) .> 0)

    LR_presence_absence = spatial_obj_exp_LR_subset_raw[LR_subset_raw_binary_mask_row, LR_subset_raw_binary_mask_col]
    LR_presence_absence_mat = LR_presence_absence.layers["binary"]

    # Filter spots based on COI and Condition
    mask = (adata.obs[:, "CompositionCluster_CC"] .∈ Ref(COI)) .& (adata.obs[:, "orig.ident"] .∈ Ref(Condition))
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
    coocur_COI_exp_specie_names = coocur_COI.var.name
    cooccur_COI_res = calculate_cooccurrence_stats(Matrix(coocur_COI_exp), coocur_COI.var.name; spp_names=true)
    println("Cooccurrence calculation ended")

    println("Summary of cooccurrence results:")
    # display(R"summary(cooccur_COI_res)")

    println("Probability table of cooccurrence results:")
    # display(R"library(cooccur); prob.table(cooccur_COI_res)")

    cooccur_res_df = cooccur_COI_res[:results]
    # Add a 'pair' column to the result DataFrame
    cooccur_res_df[!, :pair12] = string.(cooccur_res_df[!, :sp1_name], "_", cooccur_res_df[!, :sp2_name])
    cooccur_res_df[!, :pair21] = string.(cooccur_res_df[!, :sp2_name], "_", cooccur_res_df[!, :sp1_name])

    all_cooccur_pairs = Set([cooccur_res_df.pair12; cooccur_res_df.pair21])
    common_pairs = intersect(LR_pairs, all_cooccur_pairs)

    COI_enriched_LRs = DataFrame(from=String[], to=String[], correlation=Float64[], ligand_FC=Float64[], Receptor_FC=Float64[])
    pair_count = 0
    for pair in common_pairs
        pair_count += 1
        println("$pair_count / $(length(common_pairs))")

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
        pair_p = cooccur_res_df[(cooccur_res_df.pair12.==pair).|(cooccur_res_df.pair21.==pair), :p_gt][1]

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

function enriched_LRs_refactored(
    adata::AnnData, COI::Vector{String}, Condition::Vector{String},
    LR_list::Vector{String}, LR_pairs::Vector{String},
    exp_th::Real, corr_th::Real)

    println("Preparing L-R presence/absence matrix")

    # Extract the expression matrix from spatial_object
    spatial_object_exp = adata.layers["counts"]
    spatial_object_exp_norm = adata.X

    # Subset the expression matrix for the interested ligands and receptors
    spatial_obj_exp_LR_subset_raw = adata[:, in.(adata.var.name, Ref(LR_list))]

    # Binarize the expression matrix based on the expression threshold
    spatial_obj_exp_LR_subset_raw_binary = spatial_obj_exp_LR_subset_raw.layers["counts"] .> exp_th
    spatial_obj_exp_LR_subset_raw.layers["binary"] = spatial_obj_exp_LR_subset_raw_binary

    LR_subset_raw_binary_mask_col = vec(sum(spatial_obj_exp_LR_subset_raw_binary, dims=1) .> 0)
    LR_subset_raw_binary_mask_row = vec(sum(spatial_obj_exp_LR_subset_raw_binary, dims=2) .> 0)

    LR_presence_absence = spatial_obj_exp_LR_subset_raw[LR_subset_raw_binary_mask_row, LR_subset_raw_binary_mask_col]
    LR_presence_absence_mat = LR_presence_absence.layers["binary"]

    # Filter spots based on COI and Condition
    mask = (adata.obs[:, "CompositionCluster_CC"] .∈ Ref(COI)) .& (adata.obs[:, "orig.ident"] .∈ Ref(Condition))
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
    coocur_COI_exp_specie_names = coocur_COI.var.name
    cooccur_COI_res = calculate_cooccurrence_stats_refactored(Matrix(coocur_COI_exp), coocur_COI.var.name; spp_names=true)
    println("Cooccurrence calculation ended")

    println("Summary of cooccurrence results:")
    # display(R"summary(cooccur_COI_res)")

    println("Probability table of cooccurrence results:")
    # display(R"library(cooccur); prob.table(cooccur_COI_res)")

    cooccur_res_df = cooccur_COI_res.results
    # Add a 'pair' column to the result DataFrame
    cooccur_res_df[!, :pair12] = string.(cooccur_res_df.sp1_name, "_", cooccur_res_df.sp2_name)
    cooccur_res_df[!, :pair21] = string.(cooccur_res_df.sp2_name, "_", cooccur_res_df.sp1_name)

    all_cooccur_pairs = Set([cooccur_res_df.pair12; cooccur_res_df.pair21])
    common_pairs = intersect(LR_pairs, all_cooccur_pairs)

    COI_enriched_LRs = DataFrame(from=String[], to=String[], correlation=Float64[], ligand_FC=Float64[], Receptor_FC=Float64[])
    pair_count = 0
    for pair in common_pairs
        pair_count += 1
        println("$pair_count / $(length(common_pairs))")

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
        pair_p = cooccur_res_df[(cooccur_res_df.pair12.==pair).|(cooccur_res_df.pair21.==pair), :p_gt][1]

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