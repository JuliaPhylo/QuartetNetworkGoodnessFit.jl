module QuartetNetworkGoodnessFit

using DataFrames
using Distributed
using NLopt
using PhyloNetworks
using Random: seed!
using SharedArrays
using SpecialFunctions: loggamma
using Statistics: mean, median
using StatsFuns: normccdf, chisqccdf, betacdf, betaccdf

export
quarnetGoFtest!,
quarnetGoFtest,
ticr,
ticr!

include("ticr.jl")
include("quarnetGoF.jl")

# check hybrid-Lambda was downloaded (during build) and define variable hybridlambda
const depsjl_path = joinpath(@__DIR__, "..", "deps", "deps.jl")
function sethybridlambdapath(path::String) global hybridlambda = path; end
warnmsg = """Functions using hybrid-Lambda for simulations will not work. Run:
    using Pkg; Pkg.build("QuartetNetworkGoodnessFit" then restart Julia and try again.
    If the error persists, install the Hybrid-Lambda dependency manually and run:
    QuartetNetworkGoodnessFit.sethybridlambdapath("path_to_your_hybrid-Lambda_executable")
    """
if !isfile(depsjl_path)
    @warn "QuartetNetworkGoodnessFit not installed properly: no deps.jl file.\n" * warnmsg
else
    include(depsjl_path)
    try # check dependencies. `check_deps` defined in `deps.jl`
        check_deps()
    catch e
        @warn "QuartetNetworkGoodnessFit not installed properly: hybrid-Lambda executable not found.\n" * warnmsg
    end
end
end # module
