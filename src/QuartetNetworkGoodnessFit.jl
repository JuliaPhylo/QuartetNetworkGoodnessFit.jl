module QuartetNetworkGoodnessFit

using DataFrames
using Distributed
using NLopt
using PhyloCoalSimulations: simulatecoalescent
using PhyloNetworks
using Random: seed!
using SNaQ
using SharedArrays
using SpecialFunctions: loggamma
using StaticArrays
using Statistics: mean, median
using StatsFuns: normccdf, chisqccdf, betacdf, betaccdf

const PN = PhyloNetworks # for easier use of internals

export
quarnetGoFtest!,
network_expectedCF,
ticr,
ticr!

include("utils.jl")
include("quarnetconcordancefactors.jl")
include("ticr.jl")
include("quarnetGoF.jl")

end # module
