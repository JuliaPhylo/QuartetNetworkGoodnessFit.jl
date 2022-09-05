using QuartetNetworkGoodnessFit
using DataFrames, CSV
using PhyloNetworks
using Test
using Distributed

@testset "quartet network GoF tests" begin
  # each file should be its own testset, to run them all even if one has a failure
  include("test_ticr.jl")
  include("test_qnetGoF.jl")
  include("test_utils.jl")
  include("test_quarnetCFs.jl")
end
