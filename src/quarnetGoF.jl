@doc raw"""
    quarnetGoFtest!(net::HybridNetwork, df::DataFrame, optbl::Bool; quartetstat])
    quarnetGoFtest!(net::HybridNetwork, dcf::DataCF,   optbl::Bool; quartetstat])
    quarnetGoFtest(dcf::DataCF, quartetstat::Symbol)

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

# output

1. p-value of the overall goodness-of-fit test
2. test statistic (z value)
3. a dictionary of the count of p-values in each of the four category
4. a vector of outlier p-values, one for each four-taxon set
5. network (first and second versions):
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
                         quartetstat::Symbol=:LRT)
    d = readTableCF(df);
    res = quarnetGoFtest!(net, d, optbl; quartetstat=quartetstat); # order of value in results "res":
    df.p_value = res[4]             # overallpval, teststat, counts, pval, net
    return res
end

function quarnetGoFtest!(net::HybridNetwork, dcf::DataCF, optbl::Bool;
                         quartetstat::Symbol=:LRT)
    if optbl
        net = topologyMaxQPseudolik!(net,dcf);
    else
        topologyQPseudolik!(net,dcf);
    end
    res = quarnetGoFtest(dcf, quartetstat);
    return (res..., net) # (overallpval, teststat, counts, alpha, pseudolik, pval, net)
end

# @doc (@doc quarnetGoFtest!) quarnetGoFtest
function quarnetGoFtest(dcf::DataCF, quartetstat::Symbol)
    for q in dcf.quartet
        q.ngenes > 0 || error("quartet $q does not have info on number of genes")
    end
    # calculate outlier p-values
    quartetstat in [:LRT, :Qlog, :pearson] || error("$quartetstat is not a valid quartetstat option")
    pval = ( quartetstat == :LRT  ? multinom_lrt(dcf) :
            (quartetstat == :Qlog ? multinom_qlog(dcf) : multinom_pearson(dcf) ))
    # bin outlier p-values
    pcat = CategoricalArrays.cut(pval,[0, 0.01, 0.05, 0.1, 1], extend = true)
    nq = dcf.numQuartets
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
    overallpval = normccdf(teststat) # one-sided: P(Z > z)
    return (overallpval, teststat, count, pval)
end

"""
    multinom_pearson(dcf::DataCF)
    multinom_pearson(quartet::Vector{Quartet})

Calculate the vector of outlier p-values (one per four-taxon set)
using Pearson's chi-squared statistic under a multinomial distribution
for the observed concordance factors.
"""
multinom_pearson(dcf::DataCF) = multinom_pearson(dcf.quartet)
function multinom_pearson(quartet::Vector{Quartet})
    nq = length(quartet)
    pval = fill(-1.0, nq)
    for i in 1:nq
        qt = quartet[i]
        phat = qt.obsCF
        p = qt.qnet.expCF
        ngenes = qt.ngenes
        ipstat = ngenes * sum((phat .- p).^2 ./ p)
        ipval = chisqccdf(2, ipstat)
        pval[i] = ipval
    end
    return pval
end

"""
    multinom_qlog(dcf::DataCF)
    multinom_qlog(quartet::Vector{Quartet})

Calculate the vector of outlier p-values (one per four-taxon set)
using the Qlog statistic (Lorenzen, 1995),
under a multinomial distribution for the observed concordance factors.
"""
multinom_qlog(dcf::DataCF) = multinom_qlog(dcf.quartet)
function multinom_qlog(quartet::Vector{Quartet})
    nq = length(quartet)
    pval = fill(-1.0, nq)
    for i in 1:nq
        qt = quartet[i]
        phat = qt.obsCF
        p = qt.qnet.expCF
        ngenes = qt.ngenes
        mysum = 0.0
        for j in 1:3
            if phat[j] > eps(Float64)
                mysum += (phat[j]-p[j])^2/(p[j]*(log(phat[j])-log(p[j])))
            end
        end
        iqstat = 2 * ngenes * mysum
        ipval = chisqccdf(2,iqstat)
        pval[i] = ipval
    end
    return pval
end

"""
    multinom_lrt(dcf::DataCF)
    multinom_lrt(quartet::Vector{Quartet})

Calculate the vector of outlier p-values (one per four-taxon set)
using the likelihood ratio test under a multinomial distribution
for the observed concordance factors.
"""
multinom_lrt(dcf::DataCF) = multinom_lrt(dcf.quartet)
function multinom_lrt(quartet::Vector{Quartet})
    nq = length(quartet)
    pval = fill(-1.0, nq)
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
        ipval = chisqccdf(2,igstat)
        pval[i] = ipval
    end
    return pval
end

