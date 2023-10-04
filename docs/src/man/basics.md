# First Steps with ISCHIA.jl

## Setting up the Environment

Make sure that you have Julia version $\ge 1.6$ installed already. To add `ISCHIA.jl` as a package dependency, we need to install it from GitHub directly, as it hasn't been released on the Julia Registery yet. 

This can be done in either of the two ways:

```julia
julia> using Pkg

julia> Pkg.add("https://github.com/aadimator/ISCHIA.jl")
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

## Load the Data

First, let's load the required libraries.

```@repl ischia
using Muon
using RData
```

