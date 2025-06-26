# notes to maintain documentation

- built with [Documenter](https://juliadocs.github.io/Documenter.jl)
- deployed [here](https://JuliaPhylo.github.io/QuartetNetworkGoodnessFit.jl/)
  (go to `dev/` or `stable/`)
  using github and files committed to the `gh-pages` branch.

## what to update

- check Julia version in `.github/workflows/documentation.yml`
- update Documenter version in `docs/Project.toml`, and fix
  anything that the update might have broken.
- to add new source files & sections to the manual: edit
  `docs/make.jl` for sections to appear in the menu bar, and
  `docs/src/index.md` to update the main page's table of contents.

## to make a local version of the website

Open `scr/build/index.html` in your browser after doing this
from the main repo:

```shell
julia --project=docs/ -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
julia --project=docs/ --color=yes docs/make.jl
```

or from `docs/`:

```shell
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
julia --project --color=yes make.jl
```

or from `docs/` again, but interactively:
```shell
pkg> activate .
pkg> instantiate
pkg> dev ~/.julia/dev/QuartetNetworkGoodnessFit
julia> include("make.jl")
```

it will:
- test the `jldoctest` blocks of examples in the docstrings
- create or updates a `build/` directory with html files:
  the files to look at and check.
