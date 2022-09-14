```@setup gof
using QuartetNetworkGoodnessFit
ENV["COLUMNS"] = 85 # displaysize(stdout)[2] is 80 by default
```

# [goodness of fit of a candidate network](@id goodness_of_fit_1)

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
We can re-run an analysis using `net1` this time:

```@repl gof
res1 = quarnetGoFtest!(net1, qCF, true; seed=271, nsim=1000);
res1[[1,2,3]] # p-value, uncorrected z, σ
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

## empirical p-value

"Nice" networks with long branch lengths in coalescent units lead to "nice"
expected quartet concordance factors close to [100%, 0%, 0%],
when gene trees are expected to agree on a good proportion of 4-taxon sets.
On these networks and data sets with few loci, the distribution of the test
z-value may be too far from a normal approximation.
[`quarnetGoFtest!`](@ref) attempts to detect this issue. If it does, it sends
a warning suggesting to use an empirical p-value instead of using the p-value
returned by the first element of the output (e.g. `res1[1]`).

Here is example code to calculate the p-value from the empirical distribution
of simulated z-values. Note that a large number of simulations is needed for
this. We used `nsim=1000` to test network 1 above, which is a good place to start.

```@repl gof
zvalue_observed = res1[2]
zvalue_bootstrap = sort!(res1[6]) # sorted z-values simulated under the network
length(zvalue_bootstrap) # just to check we ran many simulations
using Statistics # to access "mean"
pvalue = mean(zvalue_bootstrap .>= zvalue_observed) # one-sided test: Prob(Z > z)
```

In this case we get a p-value of 1.0, which is consisten with all individual
outlier p-values being large, for *all* four-taxon sets as seen earlier:
the quartet CFs couldn't fit `net1` better.
