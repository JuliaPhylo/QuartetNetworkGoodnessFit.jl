@testset "testing GoF, multinomial distribution" begin

df = CSV.read(joinpath(dirname(Base.find_package("PhyloNetworks")),"..","examples","buckyCF.csv"))
d0 = readTableCF(df)
d = deepcopy(d0)
net3 = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");

@testset "using Pearson statistic" begin
# without optimizing branch lengths
netresult1 = quarnetGoFtest!(net3,df,false; quartetstat=:pearson, correction=:none);
@test length(netresult1) == 6
@test netresult1[2] ≈ 8.589058727506838  # z stat
@test netresult1[1] ≈ 4.384234705965304e-18   # p-value
@test df[!,:p_value] ≈ [1.2435633419824544e-11,0.0009280577186649157,0.0009280577186649157,1.2435633419824544e-11,0.007580817260552542,0.9998389684303506,0.012895416058219087,0.9438090044657973,0.9471103615208266,0.9438090044657973,0.012895416058219087,0.9471103615208266,0.9956830628718893,0.24055486863965628,0.007580817260552542]
@test net3.loglik ≈ 105.30058282649648
end

@testset "LRT statistic for 4-taxon set with small expected values" begin
# without optimizing branch lengths
netresult1 = quarnetGoFtest!(net3,d,false; correction=:none); # modified d: fills in CFs expected from net3
@test netresult1[2] ≈ 6.21966321647047 # z stat
@test netresult1[3] ≈ 1.0 # sigma
@test netresult1[1] ≈ 2.491115579898031e-10  # p-value
@test netresult1[4] ≈ [0.0024449826689709165,0.01496306673600063,0.01496306673600063,0.0024449826689709165,0.04086460431063039,0.9998541057240138,0.1901450501005025,0.8909735618259936,0.9058717147295428,0.8909735618259936,0.1901450501005025,0.9058717147295428,0.9913859984840471,0.3656465603640152,0.04086460431063039]
end

@testset "Qlog statistic for 4-taxon set with small expected values" begin
netresult1 = quarnetGoFtest(d.quartet, QuartetNetworkGoodnessFit.multinom_qlog!);
@test netresult1[1] ≈ 6.21966321647047  # z stat
@test netresult1[2] ≈ [6.188115021844278e-7,.0038172375225786225,.0038172375225786225,6.188115021844278e-7,.016885279326279496,.999846547125979,.054779352829667574,.8909734349100712,.9229400733753065,.8909734349100712,.054779352829667574,.9229400733753065,.9913859984246092,.2975827064317875,.016885279326279496]
end

@testset "check for exceptions" begin
@test_throws ErrorException ticr!(net3,d,false; quartetstat = :Qlog)
@test_throws ErrorException ticr!(net3,d,false; test = :bad);
@test_throws ErrorException quarnetGoFtest!(net3,d,false; quartetstat=:maxCF);
@test_throws ErrorException quarnetGoFtest!(net3,d,false; correction=:foo)
d.quartet[1].ngenes = -1
@test_throws ErrorException quarnetGoFtest!(net3,d,false);
end

d = deepcopy(d0)
@testset "LRT stat with dependence correction" begin
# test of expectedCF_ordered
for i in 1:15 d.quartet[i].qnet.expCF = [i+0.1, i+0.2, i+0.3]; end
expCF, taxa = QuartetNetworkGoodnessFit.expectedCF_ordered(d, net3)
taxa == ["A","B","C","D","E","O"]
expCF ≈ [13.1 13.2 13.3; 3.2 3.1 3.3; 2.2 2.1 2.3; 8.2 8.3 8.1; 10.3 10.2 10.1; 15.2 15.1 15.3; 5.2 5.1 5.3; 12.2 12.3 12.1; 9.3 9.2 9.1; 14.3 14.1 14.2; 4.1 4.3 4.2; 7.3 7.2 7.1; 1.1 1.3 1.2; 11.3 11.2 11.1; 6.1 6.2 6.3]
# test of quarnetGoFtest! with correction
netresult1 = quarnetGoFtest!(net3,d,false; seed=4321, verbose=true, nsim=1, keepfiles=true);
@test netresult1[4] ≈ [0.0024449826689709165,0.01496306673600063,0.01496306673600063,0.0024449826689709165,0.04086460431063039,0.9998541057240138,0.1901450501005025,0.8909735618259936,0.9058717147295428,0.8909735618259936,0.1901450501005025,0.9058717147295428,0.9913859984840471,0.3656465603640152,0.04086460431063039]
@test netresult1[2] ≈ 6.21966321647047 # z stat, uncorrected
@test netresult1[3] != 1.0 # sigma
hldir = filter!(f -> startswith(f, "jl_"), readdir())
@test length(hldir) == 1
@test length(readdir(hldir[1])) == 1
rm(hldir[1], recursive=true)
# similar test, but with multiple processors and other options, and interesting seed:
# for s in 1:200 netresult1 = quarnetGoFtest!(net3,d,false; seed=s, nsim=5);
#                if netresult1[3] > 3.0; @show (s, netresult1); break; end;end;
Distributed.addprocs(2)
# start with: julia -p 2 --project
# or: using Distributed; @everywhere begin; using Pkg; Pkg.activate("."); using PhyloNetworks; end
@everywhere using QuartetNetworkGoodnessFit
netresult1 = quarnetGoFtest!(net3,d,false; seed=182, nsim=5);
Distributed.rmprocs(workers())
@test netresult1[4] ≈ [0.0024449826689709165,0.01496306673600063,0.01496306673600063,0.0024449826689709165,0.04086460431063039,0.9998541057240138,0.1901450501005025,0.8909735618259936,0.9058717147295428,0.8909735618259936,0.1901450501005025,0.9058717147295428,0.9913859984840471,0.3656465603640152,0.04086460431063039]
@test netresult1[2] ≈ 6.21966321647047 # z stat, uncorrected
@test netresult1[3] ≈ 3.405362128771355 # sigma
@test netresult1[6] ≈ vcat(7.4043609719886545, repeat([-0.8885233166386386],4))
end

end
