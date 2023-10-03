"""
A mutable struct representing the output of co-occurrence analysis.

# Fields
- `results::DataFrame`: A DataFrame containing analysis results.
- `positive::Int`: The count of positive results.
- `negative::Int`: The count of negative results.
- `co_occurrences::Int`: The count of co-occurrences.
- `pairs::Int`: The total number of pairs.
- `random::Int`: The count of random results.
- `unclassifiable::Int`: The count of unclassifiable results.
- `sites::Matrix{Int}`: A matrix representing site information.
- `species::Int`: The count of species.
- `percent_sig::Float64`: The percentage of significant results.
- `spp_key::Union{Nothing, DataFrame}`: A DataFrame containing species keys (optional).
- `spp_names::Union{Vector{Int}, Vector{String}}`: A vector of species names. If names weren't provided, then it would contain numerical identifiers.
- `omitted::Union{Nothing, Int}`: The count of omitted results (optional).
- `pot_pairs::Union{Nothing, Int}`: The count of potential pairs (optional).
- `true_rand_classifier::Float64`: The true random classifier value.

This mutable struct is used to encapsulate and organize the results of co-occurrence analysis.
"""
@kwdef mutable struct CooccurOutput
    results::DataFrame
    positive::Int
    negative::Int
    co_occurrences::Int
    pairs::Int
    random::Int
    unclassifiable::Int
    sites::Matrix{Int}
    species::Int
    percent_sig::Float64
    spp_names::Union{UnitRange{Int},Vector{String}}
    spp_key::Union{Nothing,DataFrame} = nothing
    omitted::Union{Nothing,Int} = nothing
    pot_pairs::Union{Nothing,Int} = nothing
    true_rand_classifier::Float64
end



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

function create_output_dataframe(incidence, N_matrix, obs_cooccur, prob_cooccur, exp_cooccur, spp_names, spp_key, prob, only_effects)

    output = DataFrame(sp1=Integer[], sp2=Integer[], sp1_inc=Integer[], sp2_inc=Integer[],
        obs_cooccur=Integer[], prob_cooccur=Real[], exp_cooccur=Real[], p_lt=Real[], p_gt=Real[])

    @showprogress "Main Comp" for row in axes(obs_cooccur, 1)
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
    eff_standard::Bool=true
)

    if type == "spp_site"
        spp_site_mat = mat
    elseif type == "site_spp"
        spp_site_mat = transpose(mat)
    else
        error("Invalid 'type' parameter")
    end

    spp_key = spp_names ? DataFrame(num=1:size(spp_site_mat, 1), spp=specie_names) : nothing
    site_mask, N_matrix = create_site_mask_and_N_matrix(spp_site_mat, site_mask)

    tsites = size(spp_site_mat, 2)
    nspp = size(spp_site_mat, 1)
    spp_pairs = binomial(nspp, 2)

    incidence = calculate_incidence_matrix(mat, site_mask)
    prob_occur = incidence ./ N_matrix

    obs_cooccur, prob_cooccur, exp_cooccur = calculate_cooccurrence_data(spp_site_mat, site_mask, spp_pairs, prob_occur, N_matrix)

    if thresh
        n_pairs = size(prob_cooccur, 1)
        mask = exp_cooccur[:, 3] .>= 1
        prob_cooccur = prob_cooccur[mask, :]
        obs_cooccur = obs_cooccur[mask, :]
        exp_cooccur = exp_cooccur[mask, :]
        n_omitted = n_pairs - size(prob_cooccur, 1)
    end

    output = create_output_dataframe(incidence, N_matrix, obs_cooccur, prob_cooccur, exp_cooccur, spp_names, spp_key, prob, only_effects)

    true_rand = count(x -> (x.p_gt >= 0.05 && x.p_lt >= 0.05 && abs(x.obs_cooccur - x.exp_cooccur) <= (tsites * true_rand_classifier)), eachrow(output))

    cooccur_out = CooccurOutput(
        results=output,
        positive=count(x -> x.p_gt < 0.05, eachrow(output)),
        negative=count(x -> x.p_lt < 0.05, eachrow(output)),
        co_occurrences=count(x -> x.p_gt < 0.05 || x.p_lt < 0.05, eachrow(output)),
        pairs=size(output, 1),
        random=true_rand,
        unclassifiable=size(output, 1) - (true_rand + count(x -> x.p_gt < 0.05, eachrow(output)) + count(x -> x.p_lt < 0.05, eachrow(output))),
        sites=N_matrix,
        species=nspp,
        percent_sig=count(x -> x.p_gt < 0.05 || x.p_lt < 0.05, eachrow(output)) / size(output, 1) * 100,
        spp_names=1:size(spp_site_mat, 1),
        true_rand_classifier=true_rand_classifier
    )

    if spp_names
        cooccur_out.spp_key = spp_key
        cooccur_out.spp_names = specie_names
    end

    if thresh
        cooccur_out.omitted = n_omitted
        cooccur_out.pot_pairs = n_pairs
    end

    if !only_effects
        return cooccur_out
    else
        return effect_sizes(cooccur_out; standardized=eff_standard)
    end
end