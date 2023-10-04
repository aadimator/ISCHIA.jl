# Converting Seurat Objects to AnnData

The recommended way now is to use [cellgeni/sceasy](https://github.com/cellgeni/sceasy). All the other methods were giving me some kind of error. Here’s how I made it work:

## Install and Load Libraries

```r
install.packages('Seurat')
install.packages('reticulate')
devtools::install_github("cellgeni/sceasy")
```

```r
library(Seurat)
library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName') #environment where AnnData is installed
```

## Load the Seurat Object

```r
load("IBD_visium_SeuratObj_small.RData") # in my case
```

## Convert to AnnData

Make sure that `reticulate` has some python version installed, and also, you need to install `anndata` python package into that environment.

```r
sceasy::convertFormat(IBD.visium.P4, from="seurat", to="anndata", outFile='Spatial.h5ad', assay="Spatial", drop_single_values = FALSE, transfer_layers = c("counts", "data", "scale.data"))
```

## Read in Python

```python
import anndata
adata = anndata.read_h5ad("Spatial.h5ad")
```

## Read in Julia

You need to make sure that `Muon` is using the latest version from GitHub. Otherwise, it might give a `categorical` error when read in the file.

```julia
add Muon#main
```

And then, we can read in the data:

```julia
using Muon

adata = readh5ad("Spatial.h5ad")
display(first(adata.obs, 5))
display(first(adata.var, 5))
adata.obsm
```

## References

- [Convert Seurat to Scanpy h5ad - My Computational Genomic Playground (zqfang.github.io)](https://zqfang.github.io/2020-04-28-seurat2scanpy/)
- [cellgeni/sceasy: A package to help convert different single-cell data formats to each other (github.com)](https://github.com/cellgeni/sceasy)
