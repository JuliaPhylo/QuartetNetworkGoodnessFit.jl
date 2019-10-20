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

# install hybrid-Lambda and define variable hybridlambda
const depsjl_path = joinpath(@__DIR__, "..", "deps", "deps.jl")
if !isfile(depsjl_path)
    warnmsg = """
    QuartetNetworkGoodnessFit not installed properly.
    run 'using Pkg; Pkg.build("QuartetNetworkGoodnessFit")', restart Julia and try again.
    If the error persists, install the Hybrid-Lambda dependency manually.
    """
    @warn warnmsg
else
    include(depsjl_path)
    try # check dependencies. `check_deps` defined in `deps.jl`
        check_deps()
    catch e # only 1 dependency: hybrid-Lambda
        warnmsg = """
        hybrid-Lambda binary not installed (or installed properly).
        Functions that depend on gene tree simulations with hybrid-Lambda will not work.
        Install Hybrid-Lambda manually, and place the binary executable where expected...
        """
        @warn warnmsg
    end
end

end # module
