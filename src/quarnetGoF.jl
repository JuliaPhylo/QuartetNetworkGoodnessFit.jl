@doc raw"""
    quarnetGoFtest!(net::HybridNetwork, df::DataFrame, optbl::Bool; quartetstat=:LRT, correction=:simulation)
    quarnetGoFtest!(net::HybridNetwork, dcf::DataCF,   optbl::Bool; quartetstat=:LRT, correction=:simulation)
    quarnetGoFtest(dcf::DataCF, quartetstat::Symbol, correction::Symbol)
    quarnetGoFtest(outlier_pvalues::AbstractVector)

Goodness-of-fit test for the adequacy of the multispecies network coalescent,
to see if a given population or species network explains the
quartet concordance factor data adequately.
The network needs to be of level 1 at most (trees fullfil this condition),
and have branch lengths in coalescent units.
The test assumes a multinomial distribution for the observed quartet
concordance factors (CF), such that information on the number of genes for each
four-taxon set (`ngenes` field) must be present.

For each four-taxon set, an outlier p-value is obtained by comparing a
test statistic (-2log likelihood ratio by default) to a chi-square distribution
with 2 degrees of freedom (see below for other options).

The four-taxon sets are then binned into categories according to their
outlier p-values: `0-0.01`, `0.01-0.05`, `0.05-0.10`, and `0.10-1`.
Finally, a one-sided goodness-of-fit test is performed on these
binned frequencies to determine if they depart from an
expected 5% p-values being below 0.05 (versus more than 5%).
fixit: explain the correction, options `:none`, simulation sample size, etc.
fixit: fix signatures and output (all versions don't output all things)

- The first version takes a `DataFrame` where each row corresponds to a given
  four-taxon set. The data frame is modified by having an additional another column
  containing the p-values corresponding to each four-taxon set.
- The second version takes a `DataCF` object and modifies it by updating
  the expected concordance factors stored in that object.
- The last version (which all others call) assumes that the expected concordance
  factors in the `DataCF` object are correctly calculated from the test network.

# arguments

- `optbl`: when `false`, `net.loglik` is updated but branch lengths are
  taken as is. When `true`, a copy of `net` is used to conduce the test,
  with optimized branch lengths (in coalescent units) and updated loglik.
  This network is returned.

# keyword arguments

- `quartetstat`: the test statistic used to obtain an outlier p-value for
  each four-taxon set, which is then compared to a chi-squared distribution
  with 2 degrees of freedom to get a p-value.
  The default is `:LRT` for the likelihood ratio:
  ``2n_{genes} \sum_{j=1}^3 {\hat p}_j (\log{\hat p}_j - \log p_j)`` where ``p_j``
  is the quartet CF expected from the network, and ``{\hat p}_j`` is the
  quartet CF observed in the data.  
  Alternatives are `:Qlog` for the Qlog statistics (Lorenzen, 1995):
  ``2n_{genes} \sum_{j=1}^3 \frac{({\hat p}_j - p_j)^2}{p_j (\log{\hat p}_j - \log p_j)}``  
  and `:pearson` for Pearon's chi-squared statistic, which behaves poorly when
  the expected count is low (e.g. less than 5):
  ``n_{genes} \sum_{j=1}^3 \frac{({\hat p}_j - p_j)^2 }{p_j}``
- `correction`: fixit, seed::Int;
   nsim=1000::Int, verbose=true::Bool

# output

1. p-value of the overall goodness-of-fit test (corrected for dependence if requested)
2. test statistic (z value), uncorrected
3. estimated σ for the test statistic, used for the correction (1.0 if no correction)
4. a vector of outlier p-values, one for each four-taxon set
5. network (first and second versions):
   `net` with loglik field updated if `optbl` is false;
    copy of `net` with optimized branch lengths and loglik if `optbl` is true
6. in case `correction = :simulation`, vector of simulated z values
   (`nothing` if `correction = :none`).

See also: [`ticr!`](@ref).

# references

Lorenzen (1995).
A new family of goodness-of-fit statistics for discrete multivariate data.
Statistics & Probability Letters, 25(4):301-307.
doi: [10.1016/0167-7152(94)00234-8](https://doi.org/10.1016/0167-7152(94)00234-8)
"""
function quarnetGoFtest!(net::HybridNetwork,  df::DataFrame, optbl::Bool; kwargs...)
    d = readTableCF(df);
    res = quarnetGoFtest!(net, d, optbl; kwargs...);
    df[!,:p_value] .= res[4] # order in "res": overallpval, uncorrected z-value, sigma, pval, ...
    return res
end

function quarnetGoFtest!(net::HybridNetwork, dcf::DataCF, optbl::Bool;
                         quartetstat::Symbol=:LRT, correction::Symbol=:simulation,
                         seed=1234::Int, nsim=1000::Int, verbose=false::Bool, keepfiles=false::Bool)
    correction in [:simulation, :none] || error("correction ($correction) must be one of :none or :simulation")
    quartetstat in [:LRT, :Qlog, :pearson] || error("$quartetstat is not a valid quartetstat option")
    if optbl
        net = topologyMaxQPseudolik!(net,dcf);
    else
        topologyQPseudolik!(net,dcf);
    end
    outlierp_fun! = ( quartetstat ==  :LRT ? multinom_lrt! :
                     (quartetstat == :Qlog ? multinom_qlog! : multinom_pearson!))
    gof_zval, outlierpvals = quarnetGoFtest(dcf.quartet, outlierp_fun!)
    sig = 1.0 # 1 if independent: correction for dependence among 4-taxon set outlier pvalues
    sim_zvals = nothing
    if correction == :simulation
        # expCF in dcf correspond to net: assumption made below
        sig, sim_zvals = quarnetGoFtest_simulation(net, dcf, outlierp_fun!, seed, nsim, verbose, keepfiles)
    end
    overallpval = normccdf(gof_zval/sig) # one-sided: P(Z > z)
    return (overallpval, gof_zval, sig, outlierpvals, net, sim_zvals)
end

# @doc (@doc quarnetGoFtest!) quarnetGoFtest
function quarnetGoFtest(quartet::Vector{Quartet}, outlierp_fun!::Function)
   for q in quartet
        q.ngenes > 0 || error("quartet $q does not have info on number of genes")
    end
    pval = fill(-1.0, length(quartet))
    outlierp_fun!(pval, quartet) # calculate outlier p-values
    gof_zval = quarnetGoFtest(pval)
    return (gof_zval, pval)
end

function quarnetGoFtest(pval::AbstractVector)
    # uncorrected z-value: from one-sided test after binning p-values into [0,.05) and [.05,1.0]
    nsmall = count(p -> p < 0.05, pval)
    ntot = count(p -> isfinite(p), pval)
    ntot == length(pval) || @warn "$(length(pval) - ntot) outlier p-value(s) is/are Inf of NaN..."
    gof_zval = (nsmall/ntot - 0.05)/sqrt(0.0475/ntot) # 0.0475 = 0.05 * (1-0.05)
    return gof_zval
end

# WARNING: assumes dcf.quartet[i].qnet.expCF correspond to net. NOT checked.
function quarnetGoFtest_simulation(net::HybridNetwork, dcf::DataCF, outlierp_fun!::Function,
        seed::Int, nsim::Int, verbose::Bool, keepfiles::Bool)
    ngenes = ceil(Int, median(q.ngenes for q in dcf.quartet))
    # how to handle case when different #genes across quartets? need taxon set for each gene, not just ngenes.
    seed!(seed) # master seed
    hlseed = rand(1:10_000_000_000, nsim) # 1 seed for each simulation
    netHL = hybridlambdaformat(net)
    netHLquoted = "'$netHL'"
    nq = length(dcf.quartet)
    pval = fill(-1.0, nq) # to be re-used across simulations, but NOT shared between processes
    sim_zval = SharedArray{Float64}(nsim) # to be shared between processes
    # expected CFs: in dcf.quartet[i].qnet.expCF[j] for 4-taxon set i and resolution j
    # BUT countquartetsintrees might list 4-taxon sets in a different order, and might list
    # the 4 taxa within a set in a different order -> need to re-order resolutions within a set
    expCF, taxa = expectedCF_ordered(dcf, net) # 'taxa' gives the map i -> taxon[i]
    hyblamdir = mktempdir(pwd()) # temporary directory to store hybrid-Lambda output files

    @sync @distributed for irep in 1:nsim
        verbose && @info "starting replicate $irep"
        gt = joinpath(hyblamdir, "genetrees_rep$irep") # root name given to hybrid-Lambda
        gtcu = gt * "_coal_unit"                 # name of output file created by hybrid-Lambda
        hlcommand = `$hybridlambda -spcu $netHLquoted -num $ngenes -seed $(hlseed[irep]) -o $gt`
        hlout = ( verbose ? stdout : devnull )
        run(pipeline(hlcommand; stderr = hlout));
        run(`sed -i "" 's/_1//g' $gtcu`); # replaces individual names like "s5_1" into "s5"
        treelist = readMultiTopology(gtcu);
        length(treelist) == ngenes || @warn "unexpected number of gene trees, file $gtcu" # sanity check
        obsCF, t = countquartetsintrees(treelist; showprogressbar=verbose)
        irep == 1 && (taxa == t || error("different order of taxa used by countquartetsintrees"))
        outlierp_fun!(pval, (q.data for q in obsCF), expCF) # calculate outlier p-values
        sim_zval[irep] = quarnetGoFtest(pval)
        keepfiles || rm(gtcu) # delete file for replicate irep
    end
    keepfiles || rm(hyblamdir) # delete (empty) directory that had output of hybrid-Lambda
    sigma = sqrt(sum(sim_zval.^2)/nsim) # estimated sigma, assuming mean=0
    return sigma, sim_zval
end

"""
    expectedCF_ordered(dcf::DataCF, net::HybridNetwork)

Expected quartet concordance factors in `dcf`, but ordered as they would be if
output by `PhyloNetworks.countquartetsintrees`.
Output: 2-dimentional `SharedArray` (number of 4-taxon sets x 3).
`dcf.quartet[i].qnet.expCF[j]` for 4-taxon set `i` and resolution `j`
is stored in row `qi` and column `k` if `qi` is the rank of 4-taxon set `i`
(see `PhyloNetworks.quartetrank`). This rank depends on how taxa are ordered.
"""
function expectedCF_ordered(dcf::DataCF, net::HybridNetwork)
    nq = length(dcf.quartet)
    expCF = SharedArray{Float64, 2}(nq,3)
    taxa = PhyloNetworks.sort_stringasinteger!(tipLabels(net)) # same as done in countquartetsintrees
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    nCk = PhyloNetworks.nchoose1234(ntax) # matrix used to ranks 4-taxon sets
    nq == nCk[ntax+1,4] || error("dcf is assumed to contain ALL $(nCk[ntax+1,4]) four-taxon sets, but contains $nq only.")
    ptype = Vector{Int8}(undef,4) # permutation vector to sort the 4 taxa
    resperm = Vector{Int8}(undef,3) # permutation to sort the 3 resolutions accordingly
    for q in dcf.quartet
        tn = map(i->taxonnumber[i], q.taxon)
        sortperm!(ptype, tn)
        qi = PhyloNetworks.quartetrank(tn[ptype]..., nCk)
        # next: find the corresponding permutation of the 3 CFs.
        # this permutation is invariant to permuting ptype like: [a b c d] -> [c d a b] or [b a d c]
        # modify ptype to have 1 first, then only 6 permutations, corresponding to 6 permutations of 3 resolutions
        ptype .-= 0x01
        if ptype[1]==0x00
            resperm[1] = ptype[2]
            resperm[2] = ptype[3]
            resperm[3] = ptype[4]
        elseif ptype[2]==0x00 # consider ptype[[2,1,4,3]] to put 0(=1-1) first
            resperm[1] = ptype[1]
            resperm[2] = ptype[4]
            resperm[3] = ptype[3]
        elseif ptype[3]==0x00 # consider ptype[[3,4,1,2]] to put 0(=1-1) first
            resperm[1] = ptype[4]
            resperm[2] = ptype[1]
            resperm[3] = ptype[2]
        else # consider ptype[[4,3,2,1]] to put 0(=1-1) first instead of last
            resperm[1] = ptype[3]
            resperm[2] = ptype[2]
            resperm[3] = ptype[1]
        end
        expCF[qi,:] = q.qnet.expCF[resperm]
    end
    return expCF, taxa
end

"""
    multinom_pearson!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})
    multinom_pearson!(pval::AbstractVector{Float64}, obsCF, expCF::AbstractMatrix{Float64})

Calculate outlier p-values (one per four-taxon set)
using Pearson's chi-squared statistic under a multinomial distribution
for the observed concordance factors.
"""
function multinom_pearson!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})
    nq = length(quartet)
    for i in 1:nq
        qt = quartet[i]
        phat = qt.obsCF
        p = qt.qnet.expCF
        ngenes = qt.ngenes
        ipstat = ngenes * sum((phat .- p).^2 ./ p)
        pval[i] = chisqccdf(2, ipstat)
    end
    return nothing
end
function multinom_pearson!(pval::AbstractVector{Float64}, obsCF, expCF::AbstractMatrix{Float64})
    for (i,o) in enumerate(obsCF) # o = [observedCF, ngenes]
        p = expCF[i,:]
        ipstat = o[4] * sum((o[1:3] .- p).^2 ./ p)
        pval[i] = chisqccdf(2, ipstat)
    end
end


"""
    multinom_qlog!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})
    multinom_qlog!(pval::AbstractVector{Float64}, obsCF, expCF::AbstractMatrix{Float64})

Calculate outlier p-values (one per four-taxon set)
using the Qlog statistic (Lorenzen, 1995),
under a multinomial distribution for the observed concordance factors.
"""
function multinom_qlog!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})
    nq = length(quartet)
    for i in 1:nq
        qt = quartet[i]
        phat = qt.obsCF
        p = qt.qnet.expCF
        ngenes = qt.ngenes
        mysum = 0.0
        for j in 1:3
            if phat[j] > eps(Float64) && !isapprox(phat[j],p[j]; atol=1e-20)
                mysum += (phat[j]-p[j])^2/(p[j]*(log(phat[j])-log(p[j])))
            end
        end
        iqstat = 2 * ngenes * mysum
        pval[i] = chisqccdf(2,iqstat)
    end
    return nothing
end
function multinom_qlog!(pval::AbstractVector{Float64}, obsCF, expCF::AbstractMatrix{Float64})
    for (i,o) in enumerate(obsCF) # o = [observedCF, ngenes]
        p = expCF[i,:]
        mysum = 0.0
        for j in 1:3
            if o[j] > eps(Float64) && !isapprox(o[j],p[j]; atol=1e-20)
                mysum += (o[j]-p[j])^2/(p[j]*(log(o[j])-log(p[j])))
            end
        end
        iqstat = 2 * o[4] * mysum
        pval[i] = chisqccdf(2,iqstat)
    end
end

"""
    multinom_lrt!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})
    multinom_lrt!(pval::AbstractVector{Float64}, obsCF, expCF::AbstractMatrix{Float64})

Calculate outlier p-values (one per four-taxon set)
using the likelihood ratio test under a multinomial distribution
for the observed concordance factors.
"""
function multinom_lrt!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})
    nq = length(quartet)
    for i in 1:nq
        qt = quartet[i]
        phat = qt.obsCF
        p = qt.qnet.expCF
        ngenes = qt.ngenes
        mysum = 0.0
        for j in 1:3
            if phat[j] > eps(Float64)
                mysum += phat[j]*log(phat[j]/p[j])
            end
        end
        igstat = 2 * ngenes * mysum
        pval[i] = chisqccdf(2,igstat)
    end
    return nothing
end
function multinom_lrt!(pval::AbstractVector{Float64}, obsCF, expCF::AbstractMatrix{Float64})
    for (i,o) in enumerate(obsCF) # o = [observedCF, ngenes]
        p = expCF[i,:]
        mysum = 0.0
        for j in 1:3
            if o[j] > eps(Float64)
                mysum += o[j]*log(o[j]/p[j])
            end
        end
        igstat = 2 * o[4] * mysum
        pval[i] = chisqccdf(2,igstat)
    end
    return nothing
end

