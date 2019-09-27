@testset "testing GoF, multinomial distribution" begin

df = CSV.read(joinpath(dirname(Base.find_package("PhyloNetworks")),"..","examples","buckyCF.csv"))
d = readTableCF(df)
net3 = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");

@testset "using Pearson statistic" begin
# without optimizing branch lengths
netresult1 = quarnetGoFtest!(net3,df,false; quartetstat=:pearson);
@test length(netresult1) == 5
@test netresult1[2] ≈ 8.589058727506838  # z stat
@test netresult1[1] ≈ 4.384234705965304e-18   # p-value
@test netresult1[3] == Dict("[0.0, 0.01)"=>6,"[0.01, 0.05)"=>2,"[0.1, 1.0]"=>7)
@test df[!,:p_value] ≈ [1.2435633419824544e-11,0.0009280577186649157,0.0009280577186649157,1.2435633419824544e-11,0.007580817260552542,0.9998389684303506,0.012895416058219087,0.9438090044657973,0.9471103615208266,0.9438090044657973,0.012895416058219087,0.9471103615208266,0.9956830628718893,0.24055486863965628,0.007580817260552542]
@test net3.loglik ≈ 105.30058282649648
end

@testset "LRT statistic for 4-taxon set with small expected values" begin
# without optimizing branch lengths
netresult1 = quarnetGoFtest!(net3,d,false); # modified d: fills in CFs expected from net3
@test netresult1[2] ≈ 6.21966321647047 # z stat
@test netresult1[1] ≈ 2.491115579898031e-10  # p-value
@test netresult1[3] == Dict("[0.0, 0.01)"=>2,"[0.01, 0.05)"=>4,"[0.1, 1.0]"=>9)
end

@testset "Qlog statistic for 4-taxon set with small expected values" begin
netresult1 = quarnetGoFtest(d, :Qlog);
@test netresult1[2] ≈ 6.21966321647047  # z stat
@test netresult1[1] ≈ 2.491115579898031e-10   # p-value
@test netresult1[3] == Dict("[0.05, 0.1)"=>2,"[0.0, 0.01)"=>4,"[0.01, 0.05)"=>2,"[0.1, 1.0]"=>7)
end

@testset "check for exceptions" begin
@test_throws ErrorException ticr!(net3,d,false; quartetstat = :Qlog)
@test_throws ErrorException ticr!(net3,d,false; test = :bad);
@test_throws ErrorException quarnetGoFtest!(net3,d,false; quartetstat=:maxCF);
d.quartet[1].ngenes = -1
@test_throws ErrorException quarnetGoFtest!(net3,d,false);
end

end
