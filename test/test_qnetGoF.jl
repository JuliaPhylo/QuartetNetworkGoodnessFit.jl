@testset "testing GoF, multinomial distribution" begin

df = DataFrame(CSV.File(joinpath(dirname(Base.find_package("PhyloNetworks")),"..","examples","buckyCF.csv")), copycols=false)
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
@test netresult1[5].loglik ≈ 105.30058282649648
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
netresult1 = QuartetNetworkGoodnessFit.quarnetGoFtest(d.quartet, QuartetNetworkGoodnessFit.multinom_qlog!);
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
@testset "with dependence correction" begin
# test of expectedCF_ordered
for i in 1:15 d.quartet[i].qnet.expCF = [i+0.1, i+0.2, i+0.3]; end
expCF, taxa = QuartetNetworkGoodnessFit.expectedCF_ordered(d, net3)
taxa == ["A","B","C","D","E","O"]
expCF ≈ [13.1 13.2 13.3; 3.2 3.1 3.3; 2.2 2.1 2.3; 8.2 8.3 8.1; 10.3 10.2 10.1; 15.2 15.1 15.3; 5.2 5.1 5.3; 12.2 12.3 12.1; 9.3 9.2 9.1; 14.3 14.1 14.2; 4.1 4.3 4.2; 7.3 7.2 7.1; 1.1 1.3 1.2; 11.3 11.2 11.1; 6.1 6.2 6.3]
# test of quarnetGoFtest! with correction
netresult1 = quarnetGoFtest!(net3,d,false; seed=4321, verbose=true, nsim=1,
                             keepfiles=true, quartetstat=:pearson);
@test netresult1[4] ≈ [1.244e-11,.0009281,.0009281,1.244e-11,.007581,.9998,.0129,.9438,.9471,.9438,.0129,.9471,.9957,.2406,.007581] rtol=1e-4
@test netresult1[2] ≈ 8.589058727506838 # z stat, uncorrected
@test netresult1[3] != 1.0 # sigma: 0.8885233166386386 = |single simulated z|
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
netresult1 = quarnetGoFtest!(net3,d,false; seed=245, nsim=5);
@test netresult1[4] ≈ [0.0024449826689709165,0.01496306673600063,0.01496306673600063,0.0024449826689709165,0.04086460431063039,0.9998541057240138,0.1901450501005025,0.8909735618259936,0.9058717147295428,0.8909735618259936,0.1901450501005025,0.9058717147295428,0.9913859984840471,0.3656465603640152,0.04086460431063039]
@test netresult1[2] ≈ 6.21966321647047 # z stat, uncorrected
#= with some fiddling with the seed, we could get:
netresult1[3] ≈ 3.405362128771355 # sigma
netresult1[6] ≈ vcat(7.4043609719886545, repeat([-0.8885233166386386],4))
Without this fiddling, all 5 simulated z's are -0.8885233166386386 and
we get the warning about the simulated z's being far from 0
=#
@test netresult1[3] > 0.888 # sigma
@test length(netresult1[6]) == 5
@test sort(netresult1[6])[1:4] ≈ repeat([-0.8885233166386386],4)
netresult1 = (@test_logs (:warn, r"far from 0") quarnetGoFtest!(net3,d,true; seed=182, nsim=2, quartetstat=:Qlog));
# just because 2 simulated z's only, and same values bc tiny network. may break with different RNG
# note: with verbose=true, we see hybrid-lambda's warnings:
# WARNING! NOT ULTRAMETRIC!!!
# WARNING: Gene tree is not ultrametric
# ... until hybrid-lambda can accommodate non-time-consistent networks reliably
Distributed.rmprocs(workers())
@test netresult1[4] ≈ [.73,.073,.073,.73,.115,.997,.706,.997,1.,.997,.706,1.,1.,.885,.115] rtol=0.01
@test netresult1[2] ≈ -0.8885233166386386 # z stat, uncorrected

# network that caused a bug in hybrid-Lambda v0.6.2-beta, see
# https://github.com/hybridLambda/hybrid-Lambda/issues/36
net = readTopology("(((A:7.13,(B:5.98)#H18:1.15::0.79):0.1,C:7.23):0.07,((D:0.0)#H19:6.2::0.89,(E:5.64,(O:0.0,#H19:0.0::0.11):5.64):0.56):1.1,#H18:1.32::0.21);")
@test_logs quarnetGoFtest!(net,d,false; seed=419, nsim=5);

end

end
