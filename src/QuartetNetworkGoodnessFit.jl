module QuartetNetworkGoodnessFit

using CategoricalArrays: cut
using DataFrames
using NLopt
using PhyloNetworks
using SpecialFunctions: loggamma
using Statistics: mean
using StatsFuns: normccdf, chisqccdf, betacdf, betaccdf
using StatsBase: countmap

import PhyloNetworks: ticr
import PhyloNetworks: ticr!
# fixit: delete from PhyloNetworks:
# - the ticr and ticr! functions: src/ticr.jl
# - mention of ticr in PhyloNetwork.jl and in docs/src/lib/public.md
# - test_ticr.jl
# - dependencies? like: "using SpecialFunctions: loggamma, gamma"
# then delete "import PhyloNetworks: ticr, ticr!" from this package

export
quarnetGoFtest!,
quarnetGoFtest,
ticr,
ticr!

include("ticr.jl")
include("quarnetGoF.jl")

end # module
