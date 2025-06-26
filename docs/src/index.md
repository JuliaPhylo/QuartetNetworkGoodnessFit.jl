# Quarnet GoF --or simply QGoF

QuartetNetworkGoodnessFit.jl is a Julia
[package](https://github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl)
for phylogenetic networks analyses using four-taxon subsets.
It includes tools
to measure the goodness of fit
of a candidate network to data on subsets of 4 tips.
It depends on the [PhyloNetworks](https://github.com/JuliaPhylo/PhyloNetworks.jl)
package.

For a tutorial, see the manual:

```@contents
Pages = [
    "man/gof.md",
    "man/simulate.md",
    "man/expected_qCFs.md",
]
Depth = 1
```

References:

- Ruoyi Cai & Cécile Ané (2021).
  Assessing the fit of the multi-species network coalescent to multi-locus data.
  [Bioinformatics](https://doi.org/10.1093/bioinformatics/btaa863),
  37(5):634-641.
- Noah W. M. Stenz, Bret Larget, David A. Baum, and Cécile Ané (2015).
  Exploring tree-like and non-tree-like patterns using genome sequences:
  An example using the inbreeding plant species *Arabidopsis thaliana* (L.) Heynh. [Systematic Biology](https://doi.org/10.1093/sysbio/syv039), 64(5):809-823.
- [Addendum](http://www.stat.wisc.edu/~ane/publis/2015Stenz_TICR_addendum.pdf)
  describing a modification to the model in the original TICR test.
- for the simulation software and the coalescent model with correlated
  inheritance (e.g. from locus-specific gene flow probabilities):\
  John Fogg, Elizabeth S. Allman, and Cécile Ané (2023).
  PhyloCoalSimulations: A simulator for network multispecies coalescent models,
  including a new extension for the inheritance of gene flow.
  [Systematic Biology](https://doi.org/10.1093/sysbio/syad030),
  72(5):1171–1179.
- for the algorithm to get quartet concordance factors expected from a network:\
  C. Ané, J. Fogg, E.S. Allman, H. Baños and J.A. Rhodes (2024).
  Anomalous networks under the multispecies coalescent: theory and prevalence.
  [J. Math. Biol.](https://doi.org/10.1007/s00285-024-02050-7) 88:29.

## functions

```@index
Pages = ["lib/public.md", "lib/internals.md"]
Order = [:function]
```
