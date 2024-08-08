using ISCHIA
using Documenter

DocMeta.setdocmeta!(ISCHIA, :DocTestSetup, :(using ISCHIA); recursive=true)

makedocs(;
    modules=[ISCHIA],
    authors="Aadam <aadimator@gmail.com>",
    repo="https://github.com/aadimator/ISCHIA.jl/blob/{commit}{path}#{line}",
    sitename="ISCHIA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aadimator.github.io/ISCHIA.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Basic Usage" => "man/basics.md",
        "Convert Seurat Objects" => "man/seurat_conversion.md",
        "API" => Any[
            "Types" => "library/types.md",
            "Functions" => "library/functions.md",
            hide("Internals" => "library/internals.md"),
        ]
    ],
)

deploydocs(;
    repo="github.com/aadimator/ISCHIA.jl",
    devbranch="main",
)
