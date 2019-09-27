using Documenter
using QuartetNetworkGoodnessFit
DocMeta.setdocmeta!(QuartetNetworkGoodnessFit, :DocTestSetup, :(using QuartetNetworkGoodnessFit); recursive=true)

makedocs(
    sitename = "QuartetNetwork GoF.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"), # easier local build
    modules = [QuartetNetworkGoodnessFit],
    pages = [
        "home" => "index.md",
        "library" => [
            "public" => "lib/public.md",
            "internal" => "lib/internal.md",
        ]
    ],
)

deploydocs(
    repo = "github.com/cecileane/QuartetNetworkGoodnessFit.jl.git"
)
