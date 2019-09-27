using QuartetNetworkGoodnessFit
using CSV
using PhyloNetworks
using Test

@testset "quartet network GoF tests" begin
  # each file should be its own testset, to run them all even if one has a failure
  include("test_ticr.jl")
  include("test_qnetGoF.jl")
end
