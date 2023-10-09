using CSV
using Muon
using RData
using Revise
using ISCHIA
using DataFrames
using Combinatorics

datapath = "data/pdac_spatial.h5ad"

adata = readh5ad(datapath)
lr_network = load("data/lr_network.rds")

gene_names = adata.var.name

# Create LR_Pairs column
lr_network[!, :LR_Pairs] = string.(lr_network.from, "_", lr_network.to);
lr_network = lr_network[:, [:from, :to, :LR_Pairs]]

# Filter lr_network based on conditions
from_filter = in.(lr_network[:, :from], Ref(gene_names))
to_filter = in.(lr_network[:, :to], Ref(gene_names))
all_LR_network = lr_network[from_filter .& to_filter, :];

# Extract unique genes and common genes
all_LR_genes = unique(vcat(all_LR_network[:, :from], all_LR_network[:, :to]))
all_LR_genes_comm = intersect(all_LR_genes, collect(gene_names));

# Create LR.pairs and LR.pairs.AllCombos
LR_pairs = all_LR_network[:, :LR_Pairs]
all_combos = [join(combo, "_") for combo in combinations(all_LR_genes_comm, 2)];

spatial_object = adata
spatial_object.var_names = spatial_object.var.name
Condition = unique(spatial_object.obs[!, "orig.ident"])
LR_list = all_LR_genes_comm
LR_pairs = LR_pairs
exp_th = 1
corr_th = 0.2;

cc_list = ["CC3", "CC7"]

for cc in cc_list
    println("Running for $cc")
    lr_result = find_enriched_LR_pairs(spatial_object, [cc], Condition, LR_list, LR_pairs, exp_th, corr_th, cc_column="cc_ischia_10");
    CSV.write("outputs/$(cc)_lr_enrichment.csv", unique(lr_result["enriched_LRs"], :correlation))
    CSV.write("outputs/$(cc)_cooccurr_mat.csv", lr_result["cooccurrence_table"].results)
end