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

export calculate_cooccurrence_stats,
    find_enriched_LR_pairs,
    find_differentially_cooccurring_LR_pairs,
    CooccurOutput,
    summarize_cooccur

include("cooccur.jl")
include("lr_enrichment.jl")
include("utils.jl")


end # module ISCHIA
