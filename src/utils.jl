"""
    effect_sizes(cooccur; standardized=true)

Calculate effect sizes for co-occurrence analysis results.

# Arguments
- `cooccur`: A co-occurrence analysis object.
- `standardized`: A boolean indicating whether to calculate standardized effect sizes.

# Returns
Effect sizes as a data frame.

# Example
```julia
effect_sizes(cooccur_object)
"""
function effect_sizes(cooccur; standardized=true)
    cooccur_results = cooccur[:results]
    if "sp1_name" in names(cooccur_results)
        species1 = string.(cooccur_results.sp1_name)
        species2 = string.(cooccur_results.sp2_name)
    else
        species1 = cooccurrence_results.sp1
        species2 = cooccurrence_results.sp2
    end

    if standardized
        site_matrix = cooccur[:sites]
        raw_species1 = cooccur_results.sp1
        raw_species2 = cooccur_results.sp2
        raw_observed = cooccur_results.obs_cooccur
        raw_expected = cooccur_results.exp_cooccur

        effect_sizes_df = DataFrame(
            species1=String[],
            species2=String[],
            effects=Float64[]
        )

        for i in eachindex(species1)
            effect = (raw_observed[i] - raw_expected[i]) / site_matrix[raw_species1[i], raw_species2[i]]
            push!(effect_sizes_df, [species1[i], species2[i], effect])
        end

    else
        effect_sizes_df = DataFrame(
            species1=species1,
            species2=species2,
            effects=cooccur_results.obs_cooccur - cooccur_results.exp_cooccur
        )
    end

    return effect_sizes_df
end