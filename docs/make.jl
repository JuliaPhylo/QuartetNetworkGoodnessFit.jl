using Documenter, DocumenterMarkdown
using PhyloNetworks
using QuartetNetworkGoodnessFit
DocMeta.setdocmeta!(QuartetNetworkGoodnessFit, :DocTestSetup, :(using QuartetNetworkGoodnessFit); recursive=true)

makedocs(
    sitename = "QGoF",
    authors = "Cécile Ané and Ruoyi Cai",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"), # easier local build
    modules = [QuartetNetworkGoodnessFit],
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
    repo = "github.com/cecileane/QuartetNetworkGoodnessFit.jl.git",
    push_preview = true,
)
