"""
    ticr!(net::HybridNetwork, df::DataFrame, optbl::Bool; quartetstat, test)
    ticr!(net::HybridNetwork, dcf::DataCF,   optbl::Bool; quartetstat, test)
    ticr(dcf::DataCF, quartetstat::Symbol, test::Symbol)

Goodness-of-fit test for the adequacy of the multispecies network coalescent,
to see if a given population or species network explains the
quartet concordance factor data adequately. See Stenz et al 2015
and [addendum](http://www.stat.wisc.edu/~ane/publis/2015Stenz_TICR_addendum.pdf)
for the method on trees, from which the acronym TICR was derived:
Tree Incongruence Checking with R.

The tree / network needs to have branch lengths in coalescent units,
and must be of level 1. It must be fully resolved, such that the "major" quartet is
clearly defined, even though branch lengths can be 0.

The Dirichlet distribution is used to fit the distribution of the observed
quartet concordance factors (CF), with a concentration parameter estimated
from the data.

The four-taxon sets are then binned into categories according to their
outlier p-values: `0-0.01`, `0.01-0.05`, `0.05-0.10`, and `0.10-1`.
Finally, a one-sided goodness-of-fit test is performed on these
binned frequencies to determine if they depart from an
expected 5% p-values being below 0.05 (versus more than 5%),
using the default `test = :onesided`.
Alternatively, the option `test = :goodness` carries out the overall
goodness-of-fit chi-square test proposed by Stenz et al. (2015),
to test for any kind of departure from the expected frequencies
(0.01, 0.04, 0.05, 0.90) across all 4 bins.

- The first version takes a `DataFrame` where each row corresponds to a given
  four-taxon set. The data frame is modified by having an additional another column
  containing the outlier p-values corresponding to each four-taxon set.
- The second version takes a `DataCF` object and modifies it by updating
  the expected concordance factors stored in that object.
- The last version (which all others call) assumes that the expected concordance
  factors in the `DataCF` object are correctly calculated from the test network.

# arguments

- `optbl`: when `false`, the loglik field of `net` is updated but branch lengths are
  taken as is. When `true`, a copy of `net` is used to conduce the test,
  with updated branch lengths (in coalescent units) and updated loglik.
  This network is returned.
- `quartetstat = :maxCF`: test statistic used to obtain an outlier
  p-value for each four-taxon set. By default, it is the absolute difference
  between the observed CF and expected CF of the major resolution (which has
  the largest CF) if `quartetstat=:maxCF`, as used in Stenz et al. (2015).
  The other option is `quartetstat=:minpval`, in which case the absolute
  difference between the observed CF and expected CF is calculated for all
  3 resolutions, leading to 3 (non-independent) p-values. The outlier p-value
  for a given four-taxon set is taken as the minimum of these 3 p-values.
  This option can detect all kinds of departure from the network model
  for a given four-taxon set, but is not recommended because it is very liberal.
- `test = :onesided`: the overall goodness-of-fit test performed on
  binned frequencies of quartet outlier p-values; see above.

# output

1. p-value of the overall goodness-of-fit test
2. test statistic (z value)
3. a dictionary of the count of p-values in each of the four category
4. a vector of two test parameters for the Dirichlet distribution:
   - value of the concentration parameter α
   - value of the pseudo likelihood (optimized at α)
5. a vector of outlier p-values, one for each four-taxon set
6. network (first and second versions):
   `net` with loglik field updated if `optbl` is false;
    copy of `net` with optimized branch lengths and loglik if `optbl` is true

# references

NWM Stenz, B Larget, DA Baum and C Ané (2015).
Exploring tree-like and non-tree-like patterns using genome sequences:
An example using the inbreeding plant species *Arabidopsis thaliana* (L.) Heynh.
Systematic Biology, 64(5):809-823.
doi: [10.1093/sysbio/syv039](https://doi.org/10.1093/sysbio/syv039)
"""
function ticr!(net::HybridNetwork, df::DataFrame, optbl::Bool;
               quartetstat::Symbol=:maxCF, test::Symbol=:onesided)
    d = readTableCF(df);
    res = ticr!(net, d, optbl; quartetstat=quartetstat, test=test); # order of value in results "res":
    df.p_value = res[5]             # overallpval, teststat, counts, alpha_pseudolik, pval, net
    return res
end

function ticr!(net::HybridNetwork, dcf::DataCF, optbl::Bool; quartetstat::Symbol=:maxCF, test::Symbol=:onesided)
    if optbl
        net = topologyMaxQPseudolik!(net,dcf);
    else
        topologyQPseudolik!(net,dcf);
    end
    res = ticr(dcf, quartetstat, test);
    return (res..., net) # (overallpval, teststat, counts, alpha_pseudolik, pval, net)
end

# @doc (@doc ticr!) ticr
function ticr(dcf::DataCF, quartetstat::Symbol, test::Symbol)
    test in [:onesided, :goodness] || error("no valid test method specified for finding the overall p value")
    # calculate outlier p-values
    quartetstat in [:maxCF, :minpval] || error("'quartetstat' must be either need :maxCF or :minpval")
    res = ( quartetstat == :maxCF ? dirichlet_max(dcf) : dirichlet_min(dcf) )
    pval = res[1];
    testpram = [res[2],res[3]];
    # bin the outlier p-values
    counts = pval_4categorycounts(pval)
    nq = sum(values(counts)) # should be number of quartets... unless some had a NaN or Inf outlier pvalue
    e = [0.01,0.04,0.05,0.90] * nq # expected counts
    c = [counts[i] for i in ["[0.0, 0.01)","[0.01, 0.05)","[0.05, 0.1)","[0.1, 1.0]"]]
    # observed counts, but as a vector (not dictionary) and in correct order
    if test == :onesided # reduce to [0,.05) and [.05,1.0] and do a one-sided test
        teststat = ((c[1]+c[2])/nq - 0.05)/sqrt(0.0475/nq) # z-value. 0.0475 = 0.05 * (1-0.05)
        overallpval = normccdf(teststat) # one-sided: P(Z > z)
    elseif test == :goodness # as in Stenz et al (2015) and in phylolm, chi-square test statistics:
        teststat = sum((c.-e).^2 ./ e)
        overallpval = chisqccdf(3,teststat)
    end
    return (overallpval, teststat, counts, testpram, pval)
end

function pval_4categorycounts(pval::AbstractArray)
    counts = Dict( # bin p-values in 4 categories, focus on smaller p-values
        "[0.0, 0.01)" => count(p -> p < 0.01, pval),
        "[0.01, 0.05)"=> count(p -> 0.01 <= p < 0.05, pval),
        "[0.05, 0.1)" => count(p -> 0.05 <= p < 0.1, pval),
        "[0.1, 1.0]"  => count(p -> 0.1 <= p, pval))
    return counts
end


"""
    dirichlet_max(dcf::DataCF)

Calculate outlier p-values, one for each four-taxon set, using the maximum
concordance factor under a Dirichlet distribution.
Used by [`ticr!`](@ref).

# output
- vector of outlier p-values, one for each 4-taxon set
- value of the concentration parameter α
- value of the pseudo likelihood (optimized at α)
"""
function dirichlet_max(dcf::DataCF)
    res_alpha = ticr_optimalpha(dcf)
    alpha = res_alpha[2][1]
    pseudolik = res_alpha[1]
    nq = dcf.numQuartets
    pval = fill(-1.0,nq)
    for i in 1:nq
        phat = dcf.quartet[i].obsCF
        p = dcf.quartet[i].qnet.expCF
        p_max, max_idx = findmax(p)
        p_max_hat = getindex(phat,max_idx)
        p_sort = sort(dcf.quartet[i].qnet.expCF)
        # abs(p_max-p_sort[end-1]) > 1e-6 || @warn "Check the network for major quartet"
        d = abs(p_max_hat - p_max)
        temp = [1-(1-p_max)*alpha/2, 0.0]
        shapeadd = maximum(temp)
        ipval = betacdf(alpha*p_max+shapeadd, alpha*(1-p_max)+2*shapeadd, p_max-d)+
        betaccdf(alpha*p_max+shapeadd, alpha*(1-p_max)+2*shapeadd, p_max+d)
        pval[i] = ipval
    end
    return(pval, alpha, pseudolik)
end

"""
    dirichlet_min(dcf::DataCF)

First calculate outlier p-values using each of the three concordance factors
for each quartet under a Dirichlet distribution, then take the smallest
p-value among the three, as the outlier p-value for each four-taxon set.

# output
- a vector of outlier p-values, one for each quartet
- value of the concentration parameter α
- value of the pseudo likelihood (optimized at α)
"""
function dirichlet_min(dcf::DataCF)
    res_alpha = ticr_optimalpha(dcf)
    alpha = res_alpha[2][1]
    pseudolik = res_alpha[1]
    nq = dcf.numQuartets
    pval = fill(-1.0,nq)
    for i in 1:nq
        phat = dcf.quartet[i].obsCF
        p = dcf.quartet[i].qnet.expCF
        subpval = fill(-1.0,3)
        for j in 1:3 # we now calculate a p-value for every CF
            p_max = p[j]
            p_max_hat = phat[j]
            d = abs(p_max_hat - p_max)
            temp = [1-(1-p_max)*alpha/2, 0.0]
            shapeadd = maximum(temp)
            jpval = betacdf(alpha*p_max+shapeadd, alpha*(1-p_max)+2*shapeadd, p_max-d)+
            betaccdf(alpha*p_max+shapeadd, alpha*(1-p_max)+2*shapeadd, p_max+d)
            subpval[j] = jpval
        end
        pval[i] = minimum(subpval)
    end
    return(pval, alpha, pseudolik)
end

"""
    ticr_optimalpha(dcf::DataCF)

Find the concentration parameter α by maximizing the pseudo-log-likelihood
of observed quartet concordance factors.
The model assumes a Dirichlet distribution with mean equal to the expected
concordance factors calculated from a phylogenetic network (under ILS).
These expected CFs are assumed to be already calculated, and stored in `dcf`.

When calculating the pseudo-log-likelihood, this function checks the
observed concordance factors for any values equal to zero: they cause a problem
because the Dirichlet density is 0 at 0 (for concentration α > 1).
Those 0.0 observed CF values are re-set to the minimum of:
- the minimum of all expected concordance factors, and
- the minimum of all nonzero observed concordance factors.

# output
- maximized pseudo-loglikelihood
- value of α where the pseudo-loglikelihood is maximized
- return code of the optimization

The optimization uses NLOpt, with the `:LN_BOBYQA` method.
Optional arguments can tune the optimization differently:
`nloptmethod`, `xtol_rel` (1e-6 by default),
starting α value `x_start` (1.0 by default).
"""
function ticr_optimalpha(dcf::DataCF; x_start::Float64=1.0,
        nloptmethod::Symbol=:LN_BOBYQA, xtol_rel::Float64=1e-6)
    m = dcf.numQuartets
    logcftilde = 0.0
    minexpcf = minimum([minimum(q.qnet.expCF) for q in dcf.quartet])
    minnonzero_obscf = minimum([minimum(filter(!iszero,q.obsCF)) for q in dcf.quartet])
    minobscf = min(minexpcf,minnonzero_obscf)
    for q in dcf.quartet
        q.obsCF[findall(iszero, q.obsCF)] .= minobscf
    end
    for i in 1:m
        lobscf = log.(dcf.quartet[i].obsCF)
        p = dcf.quartet[i].qnet.expCF
        for i in 1:length(lobscf)
            logcftilde += lobscf[i]*p[i]
        end
    end
    function obj(x::Vector, grad::Vector)
        a = x[1] # concentration alpha
        su = 0.0
        logcfbar = 0.0
        for i in 1:m
            lobscf = log.(dcf.quartet[i].obsCF)
            p = dcf.quartet[i].qnet.expCF
            p_min = minimum(p)
            if a >= 1.0/p_min
                b = 0.0
            else
                b = 1.0 - a * p_min
            end
            logcfbar += mean(lobscf) * (1.0 - b)
            su += loggamma(a+3*b) - loggamma(a*p[1]+b) - loggamma(a*p[2]+b) - loggamma(a*p[3]+b)
        end
        pseudologLik = su + a*logcftilde -3*logcfbar
        return pseudologLik
    end
    opt = NLopt.Opt(nloptmethod,1)
    NLopt.lower_bounds!(opt, 1e-15)
    NLopt.upper_bounds!(opt, 1e5)
    NLopt.maxeval!(opt,1000)
    NLopt.xtol_rel!(opt,xtol_rel)
    NLopt.max_objective!(opt,obj)
    fmax, xmax, ret = NLopt.optimize(opt,[x_start])
    return fmax, xmax, ret
end
