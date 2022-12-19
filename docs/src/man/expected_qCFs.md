```@setup expcf
using PhyloNetworks, QuartetNetworkGoodnessFit
ENV["COLUMNS"] = 125 # displaysize(stdout)[2] is 80 by default
```
# expected concordance factors

The quartet concordance factors expected from a network can be calculated
exactly on a network of any level.
This can be done here using a recursive algorithm -- although it may run slowly
as its complexity was not optimized.
Below is an example, with a network on 6 taxa with 2 reticulations.

```@repl expcf
net = readTopology("(D:1,((C:1,#H25:0):0.1,((((B1:10,B2:1):1.5,#H1:0):10.8,((A1:1,A2:1):0.001)#H1:0::0.5):0.5)#H25:0::0.501):1);");
eCFs,t = network_expectedCF(net);
t # taxon list
first(eCFs, 2) # first 2 4-taxon sets, each with 3 quartet CFs. taxon numbers are indices in the taxon list above
```

To look at these quartet CFs, we define below a function to convert them to a string.
This string can then be printed to the screen (shown below with `print`) or written
to a file (using `write`, not shown).

```@repl expcf
qCFstring(qlist,taxa, sigdigits=4) = join(
  [join(taxa[q.taxonnumber],",") * ": " * string(round.(q.data, sigdigits=sigdigits))
   for q in qlist],   "\n");
print(qCFstring(eCFs,t))
```
For each 4-taxon set above, the 3 concordance factors are for the quartets listed
in this order: 12|34, 13|24, 14|23 where 1,2,3,4 refer to the order of the taxa
listed to the left.

We can also convert our list to a data frame:
```@repl expcf
using DataFrames
df = DataFrame(
    t1 = [t[q.taxonnumber[1]] for q in eCFs],
    t2 = [t[q.taxonnumber[2]] for q in eCFs],
    t3 = [t[q.taxonnumber[3]] for q in eCFs],
    t4 = [t[q.taxonnumber[4]] for q in eCFs],
    CF12_34 = [q.data[1] for q in eCFs],
    CF13_24 = [q.data[2] for q in eCFs],
    CF14_23 = [q.data[3] for q in eCFs],
)
```

If we wanted to compare with the observed frequency among 100 gene trees,
we could use data frame manipulations to join the two data frames:

```@repl expcf
using PhyloCoalSimulations
genetrees = simulatecoalescent(net, 100, 1); # 1000 genes, 1 individual / pop
df_sim = writeTableCF(countquartetsintrees(genetrees; showprogressbar=false)...);
first(df_sim, 2)
select!(df_sim, Not(:ngenes)); # delete column with number of genes
df_both = outerjoin(df, df_sim; on=[:t1,:t2,:t3,:t4], renamecols = "exact" => "sim")
```
