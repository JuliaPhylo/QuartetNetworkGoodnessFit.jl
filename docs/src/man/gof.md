```@setup gof
using QuartetNetworkGoodnessFit
ENV["COLUMNS"] = 85 # displaysize(stdout)[2] is 80 by default
```

# Goodness of fit of a candidate network

The example data set in this package is very small (5 taxa)
on purpose: to make this tutorial run quickly.
The data list all four-taxon sets, and the concordance factor (CF)
for each of the 3 topologies on these 4 taxa, that is, the
proportion of genes estimated to have each 4-taxon unrooted topology.

```@repl gof
using QuartetNetworkGoodnessFit, DataFrames, CSV
qCF = DataFrame(CSV.File(joinpath(dirname(pathof(QuartetNetworkGoodnessFit)), "..","test","example_qCF_5taxa.csv")), copycols=false);
qCF
```

## testing candidate networks

Let's say we have two candidate networks to test:
one that is just a tree (net0),
and the other that has one reticulation (net1).

```@repl gof
using PhyloNetworks
net1 = readTopology("((((D,C),((A,B))#H1),(#H1,E)),O);");
net0 = readTopology("((((D,C),(A,B)),E),O);"); # also: majorTree(net1)
```

`net0` has clades AB and CD sister to each other,
with E sister to ABCD and an outgroup O.
`net1` has `net0` as its major tree, and has a reticulation
going from E into the stem lineage to AB.

The goal is to test the goodness-of-fit of each topology.

It's best to start with a small number of simulations for the test,
to check that there are no issue. Then we can re-run with a larger number.
Note that the first use of a function includes time to compile the function,
so it takes longer than the second use.
Below, the option `true` means that we want to optimize the branch
lengths in the network, in coalescent units, before quantifying the
goodness-of-fit.

```@repl gof
res0 = quarnetGoFtest!(net0, qCF, true; seed=201, nsim=3);
nothing # hide
```

Before re-running with more simulations, we can update the network `net0`
to the version that has optimized branch lengths:

```@repl gof
net0 = res0[5]
```
Now we re-run the test using the option `false` to not re-optimize
branch lengths. We use `nsim=200` simulations below to make
this example faster. For a real data analysis, delete the `nsim` option
to use the default instead (1000) or specify a higher value.

```@repl gof
res0 = quarnetGoFtest!(net0, qCF, false; seed=234, nsim=200);
```

In the result, the first 3 numbers are the p-value of the overall
goodness-of-fit test, the uncorrected z-value, and the
estimate of σ to correct the z-value for dependence:

```@repl gof
res0[[1,2,3]] # p-value, uncorrected z, σ
```

Here, we have evidence that the tree `net0` does not fit the data adequately.
Note that σ=1 corresponds to independent quartet outlier p-values.
Typical estimates of σ are quite larger than 1, increase
with the number of taxa, and increase with longer branch lengths
(for a given number of taxa).

The next element is the list of outlier p-values, which was also added
to the data frame that contains our data:

```@repl gof
res0[4]
qCF[:,[:t1,:t2,:t3,:t4,:p_value]]
```

The very small overall p-value indicated an excess of outlier four-taxon sets.
Here we see what these outliers are: all four-taxon sets containing
O, E and either A or B are strong outliers (very small outlier p-values).

In fact, the CF data were simulated on `net1`, in which the
AB clade received gene flow from E.
(See below if interested in the simulation).
We can re-run an analysis using `net1` this time:

```@repl gof
res1 = quarnetGoFtest!(net1, qCF, true; seed=721, nsim=200);
res1[[1,2,3]]
```

This network is found to provide an adequate absolute fit
(as expected from how the data were simulated).

Note that after optimization of branch lengths and γs
to best fit the CF data, the network was "ultrametrized" along
its major tree to assign values to missing edge lengths.
External edges are typically missing a length in coalescent units
if there was a single individual sampled per species, for example.

```@repl gof
res1[5]
```


For more options, see [`quarnetGoFtest!`](@ref), such as for
the outlier test statistic (G or likelihood ratio test by default).

## parallel computations

For larger networks, the large number of four-taxon sets causes
the test to run more slowly. The computations can be parallelized
by providing more processors to Julia. For instance,
this can be done by starting julia with the `-p` option:

```shell
julia -p 3 # 3 worker processors, 4 processors total
```

To check for progress during the test, we can
check for a new directory with a name starting with `jl_`,
such as `jl_0CVOfE` (the end of this name is randomly generated).
Then we can check the files in this directory: their names are indicative of
the simulation replicate number that is currently being processed.
In the example below, there are 3 files in my `jl_0CVOfE` directory
(because of using 3 workers), and these files show processing
replicate numbers 106, 440 and 773 (out of the 1000 replicates by default).
Each contains 200 lines (because it contains 1 gene tree
per line, since the original data had 200 genes).
It means that the simulation-based test is at about 1/3 of its way.
By the way, these files should not be modified. At the end of the
simulation-based test, these temporary files and folder are automatically
deleted.

```shell
$ wc -l jl_0CVOfE/*
     200 jl_0CVOfE/genetrees_rep106_coal_unit
     200 jl_0CVOfE/genetrees_rep440_coal_unit
     200 jl_0CVOfE/genetrees_rep773_coal_unit
```

## quartet concordance factor simulation

The data in `qCF` were simulated using
[hybrid-Lambda](https://github.com/hybridLambda/hybrid-Lambda),
which comes with `QuartetNetworkGoodnessFit`.
In the code below, hybrid-lambda is called within julia,
asking for the simulation of 200 genes (`-num 200`).

Note that hybrid-lambda assumes that the tree or network
is ultrametric, but runs even if this assumption is not met.

```julia
net1 = readTopology("((((D:0.2,C:0.2):2.4,((A:0.4,B:0.4):1.1)#H1:1.1::0.7):2.0,(#H1:0.0::0.3,E:1.5):3.1):1.0,O:5.6);");
hl = QuartetNetworkGoodnessFit.hybridlambda # path to hybrid-lambda simulator, on local machine
net1HL = hybridlambdaformat(net1) # format for the network, to use below by hybrid-lambda
run(`$hl -spcu "((((D:0.2,C:0.2)I1:2.4,((A:0.4,B:0.4)I2:1.1)H1#0.7:1.1)I3:2.0,(H1#0.7:0.0,E:1.5)I4:3.1)I5:1.0,O:5.6)I6;" -num 200 -seed 123 -o "genetrees"`)
treelist = readMultiTopology("genetrees_coal_unit")
obsCF = writeTableCF(countquartetsintrees(treelist)...)
```
