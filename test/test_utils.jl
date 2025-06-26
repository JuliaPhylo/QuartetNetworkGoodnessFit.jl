@testset "ultrametrize!" begin
n = readnewick("((((B)#H1:::0.7,D):4.0,(#H1:::0.3,E):6.2):2.0,O);")
@test @test_logs QuartetNetworkGoodnessFit.ultrametrize!(n, true)
@test writenewick(n, round=true) == "((((B:0.0)#H1:2.2::0.7,D:2.2):4.0,(#H1:0.0::0.3,E:0.0):6.2):2.0,O:8.2);"
n = readnewick("((((B)#H1:::0.7,D):4.0,(#H1:1.0::0.3,E):2.2):2.0,O);")
@test !(@test_logs (:warn, r"^could not adjust length of major edge") QuartetNetworkGoodnessFit.ultrametrize!(n, true))
@test writenewick(n, round=true) == "((((B:0.0)#H1:0.0::0.7,D:0.0):4.0,(#H1:1.0::0.3,E:1.8):2.2):2.0,O:6.0);"
n = readnewick("((((B:2.0)#H1:0.5::0.7,D):4.0,(#H1:::0.3,E):6.2):2.0,O);")
@test !(@test_logs (:warn, r"^could not adjust length of minor edge") QuartetNetworkGoodnessFit.ultrametrize!(n, true))
@test writenewick(n, round=true) == "((((B:2.0)#H1:0.5::0.7,D:2.5):4.0,(#H1:0.0::0.3,E:0.3):6.2):2.0,O:8.5);"
n = readnewick("((((B)#H1:0.6::0.7,D):4.0,(#H1:0.1::0.3,E):6.2):2.0,O);")
@test !(@test_logs (:warn,r"hybrid H") QuartetNetworkGoodnessFit.ultrametrize!(n, true))
end

@testset "reroot!" begin
rn = readnewick("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
tn = readnewick("((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7,((#H1:0.0::0.3,E:3.0):6.2,(O:11.2):2.0):4.0);")
@test QuartetNetworkGoodnessFit.reroot!(tn, rn) == 0
tn = readnewick("((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7,((#H1:0.0::0.3,E:3.0):6.2,O:13.2):4.0);")
# same as before, except for degree-2 node removed
@test QuartetNetworkGoodnessFit.reroot!(tn, rn) == 1
end
