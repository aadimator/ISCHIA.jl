using CSV
using Muon
using RData
using Revise
using ISCHIA
using DataFrames
using Combinatorics

println("Reading data")

spatial_object = readh5ad("data/06-pdac_nac-clusters-rmv_unk.h5ad")
lr_network = load("data/lr_network.rds")

println("Data read")

gene_names = collect(spatial_object.var_names)
spatial_object.var.name = gene_names

# Create LR_Pairs column
lr_network[!, :LR_Pairs] = string.(lr_network.from, "_", lr_network.to);
lr_network = lr_network[:, [:from, :to, :LR_Pairs]];

# Filter lr_network to only include genes in adata
from_filter = in.(lr_network[!, :from], Ref(gene_names))
to_filter = in.(lr_network[:, :to], Ref(gene_names))
all_LR_network = lr_network[from_filter .& to_filter, :];

# Extract unique genes and common genes
all_LR_genes = unique(vcat(all_LR_network[:, :from], all_LR_network[:, :to]))
all_LR_genes_comm = intersect(all_LR_genes, collect(gene_names));

# Create LR.pairs and LR.pairs.AllCombos
LR_pairs = all_LR_network[:, :LR_Pairs]
all_combos = [join(combo, "_") for combo in combinations(all_LR_genes_comm, 2)];

# spatial_object = adata
LR_list = all_LR_genes_comm
LR_pairs = LR_pairs
exp_th = 1
corr_th = 0.2;

cc_column = "CC_k10"
cc_list = ["CC07"]
# Condition = unique(spatial_object.obs[!, "orig.ident"])
condition_list = ["Yes", "No"]
condition_column = "neoadjuvant_chemo"

# for cc in cc_list
#     println("Running for $cc")
#     for condition in condition_list
#         println("Running for $condition")
#         lr_result = find_enriched_LR_pairs(
#             spatial_object,
#             [cc],
#             [condition],
#             LR_list,
#             LR_pairs,
#             exp_th,
#             corr_th,
#             cc_column=cc_column,
#             condition_column=condition_column
#         )

#         # Make sure output directory exists, and if not, create it
#         if !isdir("outputs/pdac_nac")
#             mkpath("outputs/pdac_nac")
#         end

#         CSV.write("outputs/pdac_nac/$(cc)_lr_enrichment_$(condition).csv", lr_result["enriched_LRs"])
#         CSV.write("outputs/pdac_nac/$(cc)_cooccurr_mat_$(condition).csv", lr_result["cooccurrence_table"].results)
#     end
# end

# println("Running for $(cc_list)")
# println("Running for $condition_list")
lr_result = find_enriched_LR_pairs(
    spatial_object,
    cc_list,
    condition_list,
    LR_list,
    LR_pairs,
    exp_th,
    corr_th,
    cc_column=cc_column,
    condition_column=condition_column
)

# Make sure output directory exists, and if not, create it
if !isdir("outputs/pdac_nac")
    mkpath("outputs/pdac_nac")
end

CSV.write("outputs/pdac_nac/$(cc_list[1])_lr_enrichment.csv", lr_result["enriched_LRs"])

println("Done")