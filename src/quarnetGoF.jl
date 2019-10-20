@doc raw"""
    quarnetGoFtest!(net::HybridNetwork, df::DataFrame, optbl::Bool; quartetstat=:LRT, correction=:simulation)
    quarnetGoFtest!(net::HybridNetwork, dcf::DataCF,   optbl::Bool; quartetstat=:LRT, correction=:simulation)
    quarnetGoFtest(dcf::DataCF, quartetstat::Symbol, correction::Symbol)

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
- `correction`: fixit

# output

1. p-value of the overall goodness-of-fit test (corrected for dependence if requested)
2. test statistic (z value), uncorrected
3. estimated σ² for the test statistic, used for the correction (1.0 if no correction)
4. a dictionary of the count of p-values in each of the four category
5. a vector of outlier p-values, one for each four-taxon set
6. network (first and second versions):
   `net` with loglik field updated if `optbl` is false;
    copy of `net` with optimized branch lengths and loglik if `optbl` is true

See also: [`ticr!`](@ref).

# references

Lorenzen (1995).
A new family of goodness-of-fit statistics for discrete multivariate data.
Statistics & Probability Letters, 25(4):301-307.
doi: [10.1016/0167-7152(94)00234-8](https://doi.org/10.1016/0167-7152(94)00234-8)
"""
function quarnetGoFtest!(net::HybridNetwork,  df::DataFrame, optbl::Bool;
                         quartetstat::Symbol=:LRT, correction::Symbol=:simulation)
    d = readTableCF(df);
    res = quarnetGoFtest!(net, d, optbl; quartetstat=quartetstat, correction=correction);
    # order in "res": overallpval, uncorrected z-value, sigma2, counts, pval, net
    df[!,:p_value] .= res[5]
    return res
end

function quarnetGoFtest!(net::HybridNetwork, dcf::DataCF, optbl::Bool;
                         quartetstat::Symbol=:LRT, correction::Symbol=:simulation)
    if optbl
        net = topologyMaxQPseudolik!(net,dcf);
    else
        topologyQPseudolik!(net,dcf);
    end
    res = quarnetGoFtest(dcf, quartetstat, correction);
    sig = 1.0 # 1 if independent: correction for dependence among 4-taxon set outlier pvalues
    if correction == :simulation
        @warn "simulation not implemented yet"
        sig = sqrt(1.0)
    else
        @info "no correction for the dependence among the 4-taxon set outlier p-values"
    end
    overallpval = normccdf(res[1]/sig) # one-sided: P(Z > z)
    return (overallpval, res..., net)  # (overallpval, uncorrected z-value, sigma2, counts, pval, net)
end

# @doc (@doc quarnetGoFtest!) quarnetGoFtest
function quarnetGoFtest(dcf::DataCF, quartetstat::Symbol, correction::Symbol)
    for q in dcf.quartet
        q.ngenes > 0 || error("quartet $q does not have info on number of genes")
    end
    pval = fill(-1.0, dcf.numQuartets)
    quarnetGoFtest!(pval, dcf.quartet, quartetstat, correction)
end
function quarnetGoFtest!(pval::AbstractVector, quartet::Vector{Quartet}, quartetstat::Symbol, correction::Symbol)
    # calculate outlier p-values
    if quartetstat == :LRT
        multinom_lrt!(pval, quartet)
    elseif quartetstat == :Qlog
        multinom_qlog!(pval, quartet)
    elseif quartetstat == :pearson
        multinom_pearson!(pval, quartet)
    else
        error("$quartetstat is not a valid quartetstat option")
    end
    # bin outlier p-values
    nq = length(quartet)
    pcat = CategoricalArrays.cut(pval,[0, 0.01, 0.05, 0.1, 1], extend = true)
    e = [0.01,0.04,0.05,0.90] * nq # expected counts
    count = countmap(pcat)       # observed counts, but some categories might be missing
    c = Float64[]
    interval = ["[0.0, 0.01)","[0.01, 0.05)","[0.05, 0.1)","[0.1, 1.0]"]
    for i in interval
        if haskey(count, i)
            push!(c,count[i])
        else
            push!(c,0.0)
        end
    end
    # one-sided test: after reducing to [0,.05) and [.05,1.0]
    teststat = ((c[1]+c[2])/nq - 0.05)/sqrt(0.0475/nq) # z-value. 0.0475 = 0.05 * (1-0.05)
    sigma2 = 1.0
    return (teststat, sigma2, count, pval)
end

"""
    multinom_pearson!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})

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

"""
    multinom_qlog!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})

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

"""
    multinom_lrt!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})

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

