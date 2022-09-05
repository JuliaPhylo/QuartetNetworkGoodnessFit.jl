@testset "quartet CFs expected from a network" begin

qCFstring(qlist,taxa) = join([join(taxa[q.taxonnumber],",") * ": " * string(round.(q.data, sigdigits=4)) for q in qlist], "\n")

@testset "tree-like quartets, some anomalous" begin
# 4 taxa only
net = readTopology("((D:1.1,#H25:0.0):2,((((C:0.5)#H1:0.1::0.9,(#H1:0,B:0.5):0.1):0.1)#H25:0.0::0.5,A:1.1):2);")
q,t = network_expectedCF(net, showprogressbar=false);
@test qCFstring(q,t) == "A,B,C,D: [0.3721, 0.3721, 0.2559]"
#= manual calculation: major CF for AD|BC = 0.25588020130713734 from
0.1*(1-exp(-0.2)) + 0.9*(1-exp(-0.1)) + (0.1*exp(-0.2) + 0.9*exp(-0.1))*(0.5/3 + 0.5*exp(-4)/3)
minor CFs: (1-0.25588020130713734)/2 = 0.37205989934643136
=#
# 6 taxa
net = readTopology("(D:1,((C:1,#H25:0.1):0.5,((((B1:10,B2:1):1.5,#H1:0):10.8,((A1:1,A2:1):0.001)#H1:0::0.5):0.5)#H25:0.2::0.6):1);")
q,t = network_expectedCF(net); # anomalous: A1, A2, {B1 or B2}, {C or D}
# plot(net, showedgelength=true, showgamma=true)
@test t == ["A1","A2","B1","B2","C","D"]
@test qCFstring(q,t) ==
"""
A1,A2,B1,B2: [0.8885, 0.05573, 0.05573]
A1,A2,B1,C: [0.1675, 0.4162, 0.4162]
A1,A2,B2,C: [0.1675, 0.4162, 0.4162]
A1,B1,B2,C: [0.03719, 0.03719, 0.9256]
A2,B1,B2,C: [0.03719, 0.03719, 0.9256]
A1,A2,B1,D: [0.1675, 0.4162, 0.4162]
A1,A2,B2,D: [0.1675, 0.4162, 0.4162]
A1,B1,B2,D: [0.03719, 0.03719, 0.9256]
A2,B1,B2,D: [0.03719, 0.03719, 0.9256]
A1,A2,C,D: [0.6928, 0.1536, 0.1536]
A1,B1,C,D: [0.795, 0.1025, 0.1025]
A2,B1,C,D: [0.795, 0.1025, 0.1025]
A1,B2,C,D: [0.795, 0.1025, 0.1025]
A2,B2,C,D: [0.795, 0.1025, 0.1025]
B1,B2,C,D: [1.0, 9.331e-7, 9.331e-7]"""
#= manual calculations
A1A2|B1B2: (1-exp(-0.001)) + exp(-0.001)*0.25*(3*(1-exp(-1.5)*2/3) + (1-exp(-12.3)*2/3)) = 0.8885456713760765
(1-0.8885456713760765)/2 = 0.05572716431196173
AiBj|CD: 0.5*(1-exp(-0.5) + 1-exp(-11.3)) + 0.5*(exp(-0.5)+exp(-11.3))*(0.6*0.6*(1-exp(-0.2)*2/3) + 0.4*0.4*(1-exp(-0.1)*2/3) + 2*0.4*0.6*exp(-0.5)/3) = 0.7949986247417531
(1-0.7949986247417531)/2 = 0.10250068762912345
=#
end

@testset "4-blob" begin
# h=1: simple 4-cycle
net = readTopology("(D:0,((A:0)#H25:0,((#H25:0,B:0):1.2,C:0):0.9):0);")
@test_throws Exception network_expectedCF(net, showprogressbar=false); # missing γs
setGamma!(net.edge[3], 0.7) # "(D:0.0,((A:0.0)#H25:0.0::0.7,((#H25:0.0::0.3,B:0.0):1.2,C:0.0):0.9):0.0);"
net.edge[1].length=-1.0
@test_throws Exception network_expectedCF(net, showprogressbar=false); # missing edge length
net.edge[1].length=0.0
q,t = network_expectedCF(net, showprogressbar=false);
@test qCFstring(q,t) == "A,B,C,D: [0.3346, 0.125, 0.5404]"
#= manual calculation:
0.3*(1-exp(-1.2)*2/3) + 0.7*exp(-0.9)/3 = 0.3346274115570327
1-0.3346274115570327-0.5403869133122741 = 0.12498567513069325
0.3*exp(-1.2)/3 + 0.7*(1-exp(-0.9)*2/3) = 0.5403869133122741
=#
# h=2: 4-cycle + 3_2-cycle within; unrooted; 2-cycle along external edge
net = readTopology("(D:0,((((C:0)#H1:0::0.9,#H1:0):0,((B:100,#H25:100):0.01)#H22:0.01::0.8):0.05,#H22:0.2):2,(A:0)#H25:100::0.7);")
q,t = network_expectedCF(net, showprogressbar=false);
@test qCFstring(q,t) == "A,B,C,D: [0.1335, 0.1288, 0.7377]"
net.edge[10].length = 2.0 # was 0.05, within the 3_2 cycle, causing an anomaly
q,t = network_expectedCF(net, showprogressbar=false);
@test qCFstring(q,t) == "A,B,C,D: [0.08703, 0.1211, 0.7919]"
end

end
