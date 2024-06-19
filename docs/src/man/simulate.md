# quartet concordance factor simulation

## using PhyloCoalSimulations

The correction for dependence across 4-taxon uses a simulation approach.
We provide more examples here to simulate gene trees and extract the
quartet concordance factors observed in these simulated gene trees.

We use [PhyloCoalSimulations](https://JuliaPhylo.github.io/PhyloCoalSimulations.jl/stable/)
for the simulation of trees under the network multispecies coalescent,
assuming that the network edges lengths are in coalescent units.

```@repl simulate
using PhyloNetworks, QuartetNetworkGoodnessFit, PhyloCoalSimulations
net1 = readTopology("((((D:0.2,C:0.2):2.4,((A:0.4,B:0.4):1.1)#H1:1.1::0.7):2.0,(#H1:0.0::0.3,E:1.5):3.1):1.0,O:5.6);");
ngenes = 30
genetrees = simulatecoalescent(net1, ngenes, 1); # 1 individual / species
obsCF, t = countquartetsintrees(genetrees; showprogressbar=true);
df = writeTableCF(obsCF, t)
```

## using Hybrid-Lambda

!!! warning
    Hybrid-Lambda is not recommended
    (see [Allman, Baños & Rhodes 2022](https://doi.org/10.1109/TCBB.2022.3177956))
    but this section may still be of interest to some and is left for completeness.

The example `qCF` data used [earlier](@ref goodness_of_fit_1) were originally simulated as
shown below using [hybrid-Lambda](https://github.com/hybridLambda/hybrid-Lambda),
which comes with version v0.3 of `QuartetNetworkGoodnessFit`
but not with later versions.

In the code below, we show how to call hybrid-lambda within julia,
asking for the simulation of 200 genes (`-num 200`).
The major issue that this code solves is that hybrid-lambda does not
use the standard extended Newick format, but we show code to generate
the format expected by hybrid-lambda, from the standard format.

Note that hybrid-lambda assumes that the tree or network
is ultrametric, but runs even if this assumption is not met.

If you required QuartetNetworkGoodnessFit v0.3, then do this:
```julia simulate
hl = QuartetNetworkGoodnessFit.hybridlambda # path to hybrid-lambda simulator, on local machine
```
otherwise, install Hybrid-Lambda separately, then define `hl` as a string that
gives the path to the Hybrid-Lambda executable on your local machine.
Below, we use Julia to convert the standard Newick network into the format
used by Hybrid-Lambda, run it, read the gene trees from the file created by
Hybrid-Lambda, and finally calculate quartet concordance factors from these
gene trees.

```julia simulate
net1 = readTopology("((((D:0.2,C:0.2):2.4,((A:0.4,B:0.4):1.1)#H1:1.1::0.7):2.0,(#H1:0.0::0.3,E:1.5):3.1):1.0,O:5.6);");
net1HL = hybridlambdaformat(net1) # format for the network, to use below by hybrid-lambda
run(`$hl -spcu "((((D:0.2,C:0.2)I1:2.4,((A:0.4,B:0.4)I2:1.1)H1#0.7:1.1)I3:2.0,(H1#0.7:0.0,E:1.5)I4:3.1)I5:1.0,O:5.6)I6;" -num 200 -seed 123 -o "genetrees"`)
treelist = readMultiTopology("genetrees_coal_unit")
obsCF = writeTableCF(countquartetsintrees(treelist)...)
```
