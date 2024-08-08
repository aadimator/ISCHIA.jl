```@meta
CurrentModule = ISCHIA
```

# ISCHIA

This package is a reimplementation of the ISCHIA LR interaction algorithm in Julia. The original implementation is in R and can be found [here](https://github.com/ati-lz/ISCHIA). The algorithm is described in the paper [Identifying Spatial Co-occurrence in Healthy and InflAmed tissues (ISCHIA)](https://www.embopress.org/doi/full/10.1038/s44320-023-00006-5). 

In the original implementation, the LR interaction part of the algorithm is very slow and can take days to run on large datasets. This package aims to speed up the LR interaction part of the algorithm by using Julia's speed and parallelism. For comparison, the original implementation took almost 3 hours to run on their toy/smaller dataset, while the same dataset took around 6 seconds to run on this package.

Once the LR interactions are calculated, the rest of the plotting and analysis can be done in R using the original package.

```@contents
Pages = ["man/basics.md",
         "man/getting_started.md"]
Depth = 2
```
