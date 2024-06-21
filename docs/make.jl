using Documenter
using PhyloNetworks
using QuartetNetworkGoodnessFit
DocMeta.setdocmeta!(QuartetNetworkGoodnessFit, :DocTestSetup, :(using QuartetNetworkGoodnessFit); recursive=true)

makedocs(
    modules = [QuartetNetworkGoodnessFit],
    sitename = "QGoF",
    authors = "Cécile Ané and Ruoyi Cai",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true", # easier local build
        canonical="https://JuliaPhylo.github.io/QuartetNetworkGoodnessFit.jl/stable/",
        edit_link="master",
        assets=String[]
    ),
    pages = [
        "home" => "index.md",
        "manual" => [
            "goodness of fit" => "man/gof.md",
            "simulating quartet concordance factors" => "man/simulate.md",
            "expected concordance factors" => "man/expected_qCFs.md"
        ],
        "library" => [
            "public" => "lib/public.md",
            "internal" => "lib/internal.md",
        ]
    ],
)

deploydocs(
    repo = "github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl.git",
    push_preview = true,
    devbranch="master",
)
