module ISCHIA

using Muon
using RData
using Rmath
# using RCall
using Random
using DataFrames
using Statistics
using Combinatorics
using ProgressMeter

export calculate_cooccurrence_stats, enriched_LRs

"""
Calculate the co-occurrence matrix N from a binary species-site matrix.

This function creates a species by species matrix of potential co-occurring sites (N) from a binary species by site matrix, where 1 represents potential occupancy, and 0 indicates species absence.

# Arguments
- `mat::Matrix{Int}`: A binary species by site matrix.

# Returns
A species by species matrix where the upper triangle contains N for each species pair.

# Examples
```julia
# Define a binary species by site matrix
# species_matrix = rand(Bool, num_species, num_sites)

# Calculate the co-occurrence matrix N
# cooccurrence_matrix = create_N_matrix(species_matrix)

"""
function calculate_cooccurrence_matrix(mat::Matrix{Int})
    num_species = size(mat, 1)
    cooccurrence_matrix = zeros(Int, num_species, num_species)

    for i in 1:num_species
        for j in (i+1):num_species
            cooccurrence_matrix[i, j] = sum(mat[i, :] .* mat[j, :])
            cooccurrence_matrix[j, i] = cooccurrence_matrix[i, j]
        end
    end

    return cooccurrence_matrix
end

function create_site_mask_and_N_matrix(mat::Matrix{Bool}, site_mask::Union{Nothing,Matrix{Int}}=nothing)
    if !isnothing(site_mask)
        if size(site_mask) == size(mat)
            N_matrix = calculate_cooccurrence_matrix(site_mask)
        else
            error("Incorrect dimensions for 'site_mask', aborting.")
        end
    else
        site_mask = ones(Int, size(mat))
        N_matrix = size(mat, 2) * ones(Int, (size(mat, 1), size(mat, 1)))
    end

    return site_mask, N_matrix
end

function calculate_incidence_matrix(mat::Matrix, site_mask)
    nspecies = size(mat, 1)
    incidence = zeros(Int, nspecies, nspecies)
    @showprogress "Calculate Incidence" for spp in 1:nspecies
        if spp < nspecies
            for spp_next in (spp+1):nspecies
                incidence[spp, spp_next] = sum(site_mask[spp, :] .* site_mask[spp_next, :] .* mat[spp, :])
                incidence[spp_next, spp] = sum(site_mask[spp, :] .* site_mask[spp_next, :] .* mat[spp_next, :])
            end
        end
    end
    return incidence
end

function calculate_cooccurrence_data(mat::Matrix, site_mask, spp_pairs, prob_occur, N_matrix)
    obs_cooccur = zeros(Int, spp_pairs, 3)
    prob_cooccur = zeros(spp_pairs, 3)
    exp_cooccur = zeros(spp_pairs, 3)

    num_species = size(mat, 1)

    row = 0
    @showprogress "Calculate Co-occurrences" for spp in 1:num_species
        if spp < num_species
            for spp_next in (spp+1):num_species
                pairs = sum(mat[spp, site_mask[spp, :].*site_mask[spp_next, :].==1] .== 1 .&
                                                                                        mat[spp_next, site_mask[spp, :].*site_mask[spp_next, :].==1] .== 1)
                row += 1
                obs_cooccur[row, 1] = spp
                obs_cooccur[row, 2] = spp_next
                obs_cooccur[row, 3] = pairs
                prob_cooccur[row, 1] = spp
                prob_cooccur[row, 2] = spp_next
                prob_cooccur[row, 3] = prob_occur[spp, spp_next] * prob_occur[spp_next, spp]
                exp_cooccur[row, 1] = spp
                exp_cooccur[row, 2] = spp_next
                exp_cooccur[row, 3] = prob_cooccur[row, 3] * N_matrix[spp, spp_next]
            end
        end
    end
    return obs_cooccur, prob_cooccur, exp_cooccur
end

function create_output_dataframe(obs_cooccur, prob_cooccur, exp_cooccur, spp_names, spp_key, prob)

    output = DataFrame(sp1=Integer[], sp2=Integer[], sp1_inc=Integer[], sp2_inc=Integer[],
        obs_cooccur=Integer[], prob_cooccur=Real[], exp_cooccur=Real[], p_lt=Real[], p_gt=Real[])

    @showprogress "Main Comp" for row in 1:size(obs_cooccur, 1)
        sp1 = obs_cooccur[row, 1]
        sp2 = obs_cooccur[row, 2]
        sp1_inc = convert(Integer, incidence[sp1, sp2])
        sp2_inc = convert(Integer, incidence[sp2, sp1])
        max_inc = max(sp1_inc, sp2_inc)
        min_inc = min(sp1_inc, sp2_inc)
        nsite = N_matrix[sp1, sp2]
        psite = nsite + 1
        prob_share_site = zeros(Float64, psite)

        if prob == "hyper"
            if !only_effects
                all_probs = phyper.(0:min_inc, min_inc, nsite - min_inc, max_inc)
                prob_share_site[1] = all_probs[1]
                for j in 2:length(all_probs)
                    prob_share_site[j] = all_probs[j] - all_probs[j-1]
                end
            else
                for j in 0:nsite
                    if (sp1_inc + sp2_inc) <= (nsite + j)
                        if j <= min_inc
                            prob_share_site[j+1] = 1
                        end
                    end
                end
            end
        end

        if prob == "comb"
            if !only_effects
                for j in 0:nsite
                    if (sp1_inc + sp2_inc) <= (nsite + j)
                        if j <= min_inc
                            prob_share_site[j+1] = calculate_conditional_probability(j, min_inc, max_inc, nsite)
                        end
                    end
                end
            else
                for j in 0:nsite
                    if (sp1_inc + sp2_inc) <= (nsite + j)
                        if j <= min_inc
                            prob_share_site[j+1] = 1
                        end
                    end
                end
            end
        end

        p_lt = 0.0
        p_gt = 0.0
        for j in 0:nsite
            if j <= obs_cooccur[row, 3]
                p_lt += prob_share_site[j+1]
            end
            if j >= obs_cooccur[row, 3]
                p_gt += prob_share_site[j+1]
            end
            if j == obs_cooccur[row, 3]
                p_exactly_obs = prob_share_site[j+1]
            end
        end

        p_lt = round(p_lt, digits=5)
        p_gt = round(p_gt, digits=5)
        p_exactly_obs = round(p_exactly_obs, digits=5)
        prob_cooccur[row, 3] = round(prob_cooccur[row, 3], digits=3)
        exp_cooccur[row, 3] = round(exp_cooccur[row, 3], digits=1)

        push!(output, [sp1, sp2, sp1_inc, sp2_inc, obs_cooccur[row, 3],
            prob_cooccur[row, 3], exp_cooccur[row, 3], p_lt, p_gt])
    end

    if spp_names
        sp1_name = leftjoin(DataFrame(order=1:length(output.sp1), sp1=output.sp1), spp_key, on=:sp1 => :num, makeunique=true)
        sp2_name = leftjoin(DataFrame(order=1:length(output.sp2), sp2=output.sp2), spp_key, on=:sp2 => :num, makeunique=true)

        output.sp1_name = sp1_name[sortperm(sp1_name.order), "spp"]
        output.sp2_name = sp2_name[sortperm(sp2_name.order), "spp"]
    end

    return output
end

function create_output_dict(output, spp_names, specie_names, thresh, n_omitted, spp_site_mat, nspecies, true_rand_classifier)
    output_dict = Dict(:results => output,
        :positive => count(x -> x.p_gt < 0.05, eachrow(output)),
        :negative => count(x -> x.p_lt < 0.05, eachrow(output)),
        :co_occurrences => count(x -> x.p_gt < 0.05 || x.p_lt < 0.05, eachrow(output)),
        :pairs => size(output, 1),
        :random => true_rand,
        :unclassifiable => size(output, 1) - (true_rand + count(x -> x.p_gt < 0.05, eachrow(output)) + count(x -> x.p_lt < 0.05, eachrow(output))),
        :sites => N_matrix,
        :species => nspecies,
        :percent_sig => count(x -> x.p_gt < 0.05 || x.p_lt < 0.05, eachrow(output)) / size(output, 1) * 100
    )

    if spp_names
        output_dict[:spp_key] = spp_key
        output_dict[:spp_names] = specie_names
    else
        output_dict[:spp_names] = 1:size(spp_site_mat, 1)
    end

    if thresh
        output_dict[:omitted] = n_omitted
        output_dict[:pot_pairs] = n_pairs
    end
    return output_dict
end

"""
Calculate co-occurrence statistics and probabilities.

# Arguments
- `mat::Matrix{Int}`: A binary species by site matrix.
- `specie_names::Vector{String}`: Names of species.
- `type::String`: Type of matrix ('spp_site' or 'site_spp').
- `thresh::Bool`: Whether to apply a threshold.
- `spp_names::Bool`: Whether to include species names.
- `true_rand_classifier::Float64`: True random classifier.
- `prob::String`: Probability calculation method ('hyper' or 'comb').
- `site_mask::Union{Nothing, Matrix{Int}}`: Matrix specifying sites.
- `only_effects::Bool`: Whether to calculate only effects.
- `eff_standard::Bool`: Whether to standardize effects.
- `eff_matrix::Bool`: Whether to calculate effect matrix.

# Returns
A dictionary containing various co-occurrence statistics and results.
"""
function calculate_cooccurrence_stats(
    mat::Matrix{Bool},
    specie_names::Vector{String};
    type::String="spp_site",
    thresh::Bool=true,
    spp_names::Bool=false,
    true_rand_classifier::Float64=0.1,
    prob::String="hyper",
    site_mask::Union{Nothing,Matrix{Int}}=nothing,
    only_effects::Bool=false,
    eff_standard::Bool=true,
    eff_matrix::Bool=false
)
    if type == "spp_site"
        spp_site_mat = mat
    elseif type == "site_spp"
        spp_site_mat = transpose(mat)
    else
        error("Invalid 'type' parameter")
    end

    if spp_names
        spp_key = DataFrame(num=1:size(spp_site_mat, 1), spp=specie_names)
    end

    if !isnothing(site_mask)
        if size(site_mask) == size(spp_site_mat)
            N_matrix = calculate_cooccurrence_matrix(site_mask)
        else
            error("Incorrect dimensions for 'site_mask', aborting.")
        end
    else
        site_mask = ones(Int, size(spp_site_mat))
        N_matrix = size(spp_site_mat, 2) * ones(Int, (size(spp_site_mat, 1), size(spp_site_mat, 1)))
    end

    tsites = size(spp_site_mat, 2)
    nspp = size(spp_site_mat, 1)
    spp_pairs = binomial(nspp, 2)

    incidence = zeros(Int, size(N_matrix))
    prob_occur = zeros(size(N_matrix))

    obs_cooccur = zeros(Int, spp_pairs, 3)
    prob_cooccur = zeros(spp_pairs, 3)
    exp_cooccur = zeros(spp_pairs, 3)

    mat_matrix = Matrix(mat)
    @showprogress "Calculate Incidence" for spp in 1:nspp
        if spp < nspp
            for spp_next in (spp+1):nspp
                incidence[spp, spp_next] = sum(site_mask[spp, :] .* site_mask[spp_next, :] .* mat_matrix[spp, :])
                incidence[spp_next, spp] = sum(site_mask[spp, :] .* site_mask[spp_next, :] .* mat_matrix[spp_next, :])
            end
        end
    end

    prob_occur .= incidence ./ N_matrix

    row = 0
    @showprogress "Calculate Co-occurrences" for spp in 1:nspp
        if spp < nspp
            for spp_next in (spp+1):nspp
                pairs = sum(mat_matrix[spp, site_mask[spp, :].*site_mask[spp_next, :].==1] .== 1 .&
                                                                                               mat_matrix[spp_next, site_mask[spp, :].*site_mask[spp_next, :].==1] .== 1)
                row += 1
                obs_cooccur[row, 1] = spp
                obs_cooccur[row, 2] = spp_next
                obs_cooccur[row, 3] = pairs
                prob_cooccur[row, 1] = spp
                prob_cooccur[row, 2] = spp_next
                prob_cooccur[row, 3] = prob_occur[spp, spp_next] * prob_occur[spp_next, spp]
                exp_cooccur[row, 1] = spp
                exp_cooccur[row, 2] = spp_next
                exp_cooccur[row, 3] = prob_cooccur[row, 3] * N_matrix[spp, spp_next]
            end
        end
    end

    if thresh
        n_pairs = size(prob_cooccur, 1)
        mask = exp_cooccur[:, 3] .>= 1
        prob_cooccur = prob_cooccur[mask, :]
        obs_cooccur = obs_cooccur[mask, :]
        exp_cooccur = exp_cooccur[mask, :]
        n_omitted = n_pairs - size(prob_cooccur, 1)
    end

    output = DataFrame(sp1=Integer[], sp2=Integer[], sp1_inc=Integer[], sp2_inc=Integer[],
        obs_cooccur=Integer[], prob_cooccur=Real[], exp_cooccur=Real[], p_lt=Real[], p_gt=Real[])

    @showprogress "Main Comp" for row in 1:size(obs_cooccur, 1)
        sp1 = obs_cooccur[row, 1]
        sp2 = obs_cooccur[row, 2]
        sp1_inc = convert(Integer, incidence[sp1, sp2])
        sp2_inc = convert(Integer, incidence[sp2, sp1])
        max_inc = max(sp1_inc, sp2_inc)
        min_inc = min(sp1_inc, sp2_inc)
        nsite = N_matrix[sp1, sp2]
        psite = nsite + 1
        prob_share_site = zeros(Float64, psite)

        if prob == "hyper"
            if !only_effects
                all_probs = phyper.(0:min_inc, min_inc, nsite - min_inc, max_inc)
                prob_share_site[1] = all_probs[1]
                for j in 2:length(all_probs)
                    prob_share_site[j] = all_probs[j] - all_probs[j-1]
                end
            else
                for j in 0:nsite
                    if (sp1_inc + sp2_inc) <= (nsite + j)
                        if j <= min_inc
                            prob_share_site[j+1] = 1
                        end
                    end
                end
            end
        end

        if prob == "comb"
            if !only_effects
                for j in 0:nsite
                    if (sp1_inc + sp2_inc) <= (nsite + j)
                        if j <= min_inc
                            prob_share_site[j+1] = calculate_conditional_probability(j, min_inc, max_inc, nsite)
                        end
                    end
                end
            else
                for j in 0:nsite
                    if (sp1_inc + sp2_inc) <= (nsite + j)
                        if j <= min_inc
                            prob_share_site[j+1] = 1
                        end
                    end
                end
            end
        end

        p_lt = 0.0
        p_gt = 0.0
        for j in 0:nsite
            if j <= obs_cooccur[row, 3]
                p_lt += prob_share_site[j+1]
            end
            if j >= obs_cooccur[row, 3]
                p_gt += prob_share_site[j+1]
            end
            if j == obs_cooccur[row, 3]
                p_exactly_obs = prob_share_site[j+1]
            end
        end

        p_lt = round(p_lt, digits=5)
        p_gt = round(p_gt, digits=5)
        p_exactly_obs = round(p_exactly_obs, digits=5)
        prob_cooccur[row, 3] = round(prob_cooccur[row, 3], digits=3)
        exp_cooccur[row, 3] = round(exp_cooccur[row, 3], digits=1)

        push!(output, [sp1, sp2, sp1_inc, sp2_inc, obs_cooccur[row, 3],
            prob_cooccur[row, 3], exp_cooccur[row, 3], p_lt, p_gt])
    end

    if spp_names
        sp1_name = leftjoin(DataFrame(order=1:length(output.sp1), sp1=output.sp1), spp_key, on=:sp1 => :num, makeunique=true)
        sp2_name = leftjoin(DataFrame(order=1:length(output.sp2), sp2=output.sp2), spp_key, on=:sp2 => :num, makeunique=true)

        output.sp1_name = sp1_name[sortperm(sp1_name.order), "spp"]
        output.sp2_name = sp2_name[sortperm(sp2_name.order), "spp"]
    end

    true_rand = count(x -> (x.p_gt >= 0.05 && x.p_lt >= 0.05 && abs(x.obs_cooccur - x.exp_cooccur) <= (tsites * true_rand_classifier)), eachrow(output))

    output_dict = Dict(:results => output,
        :positive => count(x -> x.p_gt < 0.05, eachrow(output)),
        :negative => count(x -> x.p_lt < 0.05, eachrow(output)),
        :co_occurrences => count(x -> x.p_gt < 0.05 || x.p_lt < 0.05, eachrow(output)),
        :pairs => size(output, 1),
        :random => true_rand,
        :unclassifiable => size(output, 1) - (true_rand + count(x -> x.p_gt < 0.05, eachrow(output)) + count(x -> x.p_lt < 0.05, eachrow(output))),
        :sites => N_matrix,
        :species => nspp,
        :percent_sig => count(x -> x.p_gt < 0.05 || x.p_lt < 0.05, eachrow(output)) / size(output, 1) * 100,
        :true_rand_classifier => true_rand_classifier)

    if spp_names
        output_dict[:spp_key] = spp_key
        output_dict[:spp_names] = specie_names
    else
        output_dict[:spp_names] = 1:size(spp_site_mat, 1)
    end

    if thresh
        output_dict[:omitted] = n_omitted
        output_dict[:pot_pairs] = n_pairs
    end

    return output_dict
    # if !only_effects
    #     return output_dict
    # else
    #     return effect_sizes(output_dict, standardized=eff_standard, matrix=eff_matrix)
    # end
end


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



end # module ISCHIA
