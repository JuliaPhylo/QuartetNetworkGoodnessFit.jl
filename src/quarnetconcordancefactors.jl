"""
    network_expectedCF(net::HybridNetwork; showprogressbar=true,
            inheritancecorrelation=0)

Calculate the quartet concordance factors (qCF) expected from the multispecies
coalescent along network `net`, using the recursive algorithm by
[Ané et al. (2024)](https://doi.org/10.1007/s00285-024-02050-7).

Output: `(q,t)` where `t` is a list of taxa,
and `q` is a list of 4-taxon set objects of type `PhyloNetworks.QuartetT{datatype}`.
In each element of `q`, `taxonnumber` gives the indices in `taxa`
of the 4 taxa of interest; and `data` contains the 3 concordance factors, for the
3 unrooted topologies in the following order:
`t1,t2|t3,t4`, `t1,t3|t2,t4` and `t1,t4|t2,t3`.
This output is similar to that of `PhyloNetworks.countquartetsintrees` when
1 individual = 1 taxon, with 4-taxon sets listed in the same order
(same output `t`, then same order of 4-taxon sets in `q`).

Assumption: the network should have **edge lengths in coalescent units**.

By default, lineages at a hybrid node come from a parent (chosen according
to inheritance probabilities γ) *independently* across lineages.
With option `inheritancecorrelation > 0`, lineages have positive dependence,
e.g. to model locus-specific inheritance probabilities, randomly drawn from a
Beta distribution with mean γ across all loci
(see [Fogg, Allman & Ané 2023](https://doi.org/10.1093/sysbio/syad030)).
If `inheritancecorrelation` is
set to 1, then all lineages at a given locus inherit from the same
(randomly sampled) parent. More generally, the lineages' parents
are distributed according to a Dirichlet process with base distribution determined
by the γ values, and with concentration parameter α = (1-r)/r, that is, r = 1/(1+α),
where `r` is the input inheritance correlation.

# example

Below we use a network with two "3₂" cycles, causing some quartets to
be anomalous (see [Ané et al. 2024](https://doi.org/10.1007/s00285-024-02050-7)),
in the sense that the unrooted tree topology displayed in the network, for
example A1,A2|B1,C has a lower frequency among gene trees than another topology
not even displayed in the network.
This is due to incomplete lineage sorting interacting with gene flow.

```jldoctest
julia> using PhyloNetworks, QuartetNetworkGoodnessFit

julia> net = readnewick("(D:1,((C:1,#H25:0):0.1,((((B1:10,B2:1):1.5,#H1:0):10.8,
                ((A1:1,A2:1):0.001)#H1:0::0.5):0.5)#H25:0::0.501):1);");

julia> # using PhyloPlots; plot(net, showedgelength=true);

julia> q,t = network_expectedCF(net); # anomalous: A1, A2, {B1 or B2}, {C or D}
Calculation quartet CFs for 15 quartets...
0+---------------+100%
  ***************

julia> show(q[1].taxonnumber)
[1, 2, 3, 4]
julia> show(q[1].data)
[0.8885456713760765, 0.05572716431196174, 0.05572716431196172]

julia> for qi in q
         println(join(t[qi.taxonnumber],",") * ": " *
            string(round.(qi.data, sigdigits=3)))
       end
A1,A2,B1,B2: [0.889, 0.0557, 0.0557]
A1,A2,B1,C: [0.168, 0.416, 0.416]
A1,A2,B2,C: [0.168, 0.416, 0.416]
A1,B1,B2,C: [0.0372, 0.0372, 0.926]
A2,B1,B2,C: [0.0372, 0.0372, 0.926]
A1,A2,B1,D: [0.168, 0.416, 0.416]
A1,A2,B2,D: [0.168, 0.416, 0.416]
A1,B1,B2,D: [0.0372, 0.0372, 0.926]
A2,B1,B2,D: [0.0372, 0.0372, 0.926]
A1,A2,C,D: [0.69, 0.155, 0.155]
A1,B1,C,D: [0.793, 0.103, 0.103]
A2,B1,C,D: [0.793, 0.103, 0.103]
A1,B2,C,D: [0.793, 0.103, 0.103]
A2,B2,C,D: [0.793, 0.103, 0.103]
B1,B2,C,D: [1.0, 9.42e-7, 9.42e-7]

```
"""
function network_expectedCF(
    net::HybridNetwork;
    showprogressbar=true,
    inheritancecorrelation=0
)
    getroot(net).leaf && error("The root can't be a leaf.")
    PN.check_nonmissing_nonnegative_edgelengths(net,
        "Edge lengths are needed in coalescent units to calcualte expected CFs.")
    all(e.gamma >= 0.0 for e in net.edge) ||
        error("some γ's are missing for hybrid edges: can't calculate expected CFs.")
    inheritancecorrelation >= 0 ||
        error("the inheritance correlation should be non-negative")
    inheritancecorrelation <= 1 ||
        error("the inheritance correlation should be <= 1")
    taxa = sort!(tiplabels(net))
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    nCk = PN.nchoose1234(ntax) # matrix to rank 4-taxon sets
    qtype = MVector{3,Float64} # 3 floats: CF12_34, CF13_24, CF14_23; initialized at 0.0
    numq = nCk[ntax+1,4]
    quartet = Vector{PN.QuartetT{qtype}}(undef, numq)
    ts = [1,2,3,4]
    for qi in 1:numq
        quartet[qi] = PN.QuartetT(qi, SVector{4}(ts), MVector(0.,0.,0.))
        # next: find the 4-taxon set with the next rank,
        #       faster than using the direct mapping function
        ind = findfirst(x -> x>1, diff(ts))
        if ind === nothing ind = 4; end
        ts[ind] += 1
        for j in 1:(ind-1)
            ts[j] = j
        end
    end
    if showprogressbar
        nstars = (numq < 50 ? numq : 50)
        nquarnets_perstar = (numq/nstars)
        println("Calculation quartet CFs for $numq quartets...")
        print("0+" * "-"^nstars * "+100%\n  ")
        stars = 0
        nextstar = Integer(ceil(nquarnets_perstar))
    end

    # SNaQ takes inheritance correlation in the form of α, not ρ
    α::Float64 = inheritancecorrelation == 0.0 ? Inf : (1 - inheritancecorrelation) / inheritancecorrelation
    
    for qi in 1:numq
        quartet[qi].data .= SNaQ.compute_eCF_4taxa(net, taxa[quartet[qi].taxonnumber], α)
        
        if showprogressbar && qi >= nextstar
            print("*")
            stars += 1
            nextstar = Integer(ceil((stars+1) * nquarnets_perstar))
        end
    end
    showprogressbar && print("\n")
    return quartet, taxa
end

