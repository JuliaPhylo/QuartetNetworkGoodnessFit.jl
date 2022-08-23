# QuartetNetworkGoodnessFit.jl


Julia
[package](https://github.com/cecileane/QuartetNetworkGoodnessFit.jl)
to measure the goodness of fit
of a phylogenetic network to data on subsets of 4 tips.
It depends on the [PhyloNetworks](https://github.com/crsl4/PhyloNetworks.jl)
package.

For a tutorial, see the manual:

```@contents
Pages = [
    "man/gof.md",
    "man/simulate.md",
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

## functions

```@index
Pages = ["lib/public.md", "lib/internals.md"]
Order = [:function]
```
