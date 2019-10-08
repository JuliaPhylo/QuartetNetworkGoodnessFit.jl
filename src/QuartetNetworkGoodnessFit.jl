module QuartetNetworkGoodnessFit

using CategoricalArrays: cut
using DataFrames
using NLopt
using PhyloNetworks
using SpecialFunctions: loggamma
using Statistics: mean
using StatsFuns: normccdf, chisqccdf, betacdf, betaccdf
using StatsBase: countmap

export
quarnetGoFtest!,
quarnetGoFtest,
ticr,
ticr!

include("ticr.jl")
include("quarnetGoF.jl")

end # module
