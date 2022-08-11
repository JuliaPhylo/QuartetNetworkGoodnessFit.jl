module QuartetNetworkGoodnessFit

using DataFrames
using Distributed
using NLopt
using PhyloCoalSimulations: simulatecoalescent
using PhyloNetworks
using Random: seed!
using SharedArrays
using SpecialFunctions: loggamma
using Statistics: mean, median
using StatsFuns: normccdf, chisqccdf, betacdf, betaccdf

export
quarnetGoFtest!,
ticr,
ticr!

include("utils.jl")
include("ticr.jl")
include("quarnetGoF.jl")

end # module
