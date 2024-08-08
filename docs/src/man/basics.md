# First Steps with ISCHIA.jl

## Setting up the Environment

Make sure that you have Julia version $\ge 1.6$ installed already. To add `ISCHIA.jl` as a package dependency, we need to install it from GitHub directly, as it hasn't been released on the Julia Registery yet. 

This can be done in either of the two ways:

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/aadimator/ISCHIA.jl")
```

or

```julia
julia> ] # ']' should be pressed

(@v1.9) pkg> add https://github.com/aadimator/ISCHIA.jl
```

This will clone the repository from the GitHub and install it in your package environment.

If you already have the GitHub repository cloned on your system and you want to install the package using that, you can do so with the following commands. Make sure that you are in the cloned project directory when running the followin commands.

```julia
julia> using Pkg

julia> Pkg.add(".") # Install current directory as a package
```

or

```julia
julia> ] # ']' should be pressed

(@v1.9) pkg> add .
```

Throughout the rest of the tutorial we will assume that you have installed the ISCHIA.jl package and have already typed `using ISCHIA` which loads the package:

```jldoctest ischia
julia> using ISCHIA
```

`ISCHIA` provides optimized implementations of several functions from the [original ISCHIA](https://github.com/ati-lz/ISCHIA) repository, which was coded in `R`. Below, we'll try to replicate their **Cooccurrence of ligand and receptors** section, using their provided small dataset. Instructions on converting and importing their `Seurat Object` into `AnnData` object has been described in section.

## Load the Libraries

First, let's load the required libraries. Install any library that isn't in your package environment already. Julia will automatically suggest the command to install missing packages.

```@example ischia
using Muon
using RData
using ISCHIA
using DataFrames
using Combinatorics
```

## Load the Data

```@example ischia
sdata = readh5ad(joinpath(pkgdir(ISCHIA), "docs", "src", "assets", "Spatial.h5ad"))
lr_network = load(joinpath(pkgdir(ISCHIA), "docs", "src", "assets", "lr_network.rds"))
size(lr_network)
```

## Filter LR Network

```@example ischia
gene_names = sdata.var.name
sdata.var_names = gene_names

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
length(all_combos)
```

## Calculate LR Enrchiment

```@example ischia
spatial_object = sdata
spatial_object.var_names = spatial_object.var.name
Condition = unique(spatial_object.obs[!, "orig.ident"])
LR_list = all_LR_genes_comm
LR_pairs = LR_pairs
exp_th = 1
corr_th = 0.2

cc4_enriched_lr = find_enriched_LR_pairs(spatial_object, ["CC4"], Condition, LR_list, LR_pairs, exp_th, corr_th);
cc7_enriched_lr = find_enriched_LR_pairs(spatial_object, ["CC7"], Condition, LR_list, LR_pairs, exp_th, corr_th);
```

## Display Calculated LR Encrichment Scores

```@example ischia
first(cc4_enriched_lr["enriched_LRs"], 5)
```

```@example ischia
first(cc7_enriched_lr["enriched_LRs"], 5)
```

## Find Differentially Occurring LR pairs

```@example ischia
diff_df47 = find_differentially_cooccurring_LR_pairs(cc7_enriched_lr, cc4_enriched_lr, 0.05, 0.1)
first(diff_df47, 10)
```
