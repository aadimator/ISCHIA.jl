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
    enriched_LRs,
    calculate_cooccurrence_stats_refactored,
    enriched_LRs_refactored,
    CooccurOutput,
    summarize_cooccur

include("cooccur.jl")
include("lr_enrichment.jl")
include("utils.jl")


end # module ISCHIA
