@testset "testing TICR, Dirichlet distribution" begin
# previously in PhyloNetworks

df = CSV.read(joinpath(dirname(Base.find_package("PhyloNetworks")),"..","examples","buckyCF.csv"));
d = readTableCF(df);

@testset "ticr! on data frame, tree" begin
net1 = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
setGamma!(net1.edge[9],0.0);
net2 = deepcopy(net1);
# without optimizing branch lengthes
result1 = ticr!(net1,df,false);
@test result1[2] ≈ 1.480872194397731   # z statistic, from R: prop.test(2,15, p=0.05, alternative="greater", correct=F)
@test result1[1] ≈ 0.06932031690660927 # p-value, from R
@test result1[3] == Dict("[0.0, 0.01)" => 2, "[0.05, 0.1)" => 2, "[0.1, 1.0]"  => 11)
@test result1[4][2] ≈ 48.152697007372566418 # pseudo log-lik obtained from R
@test result1[4][1] ≈ 10.576940922426542713 atol=1e-5 # alpha obtained from R
#with optimizing branch lengthes
result2 = ticr!(net2,d,true);
net2_1 = readTopology("(D,C,(((A,B):1.297)#H1:1.412::1.0,((#H1:0.0::0.0,E):6.877,O):1.354):9.352);");
result3 = ticr!(net2_1,d,false);
@test result3[2] ≈ 1.480872194397731   # z statistic, from R, same as above
@test result3[1] ≈ 0.06932031690660927 # p-value, from R
@test result3[3] == Dict("[0.0, 0.01)" => 2, "[0.05, 0.1)" => 2, "[0.1, 1.0]"  => 11)
@test result3[4][2] ≈ 54.241562916699216 # pseudo log-lik obtained from R
@test result3[4][1] ≈ 20.128258663235194 # alpha obtained from R
result3_1 = ticr!(net2_1,d,false; test=:goodness);
@test result3_1[2] ≈ 25.962962962962965463 # chi-squared statistic obtained from R
@test result3_1[1] ≈ 9.7092282251534852702e-06 # p-value obtained from R
@test result3_1[3] == Dict("[0.0, 0.01)" => 2, "[0.05, 0.1)" => 2, "[0.1, 1.0]"  => 11)
@test result3[4][2] ≈ 54.241562916699216 atol=1e-3 # pseudo log-lik obtained from R
@test result3[4][1] ≈ 20.128258663235194 atol=1e-4 # alpha obtained from R
end

@testset "ticr! on data frame, network, p value from maximum observed values, dirichlet distribution" begin
net3 = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
# without optimizing branch lengths
netresult1 = ticr!(net3,d,false);
@test netresult1[2] ≈ -0.8885233166386386  # z stat, from R: prop.test(0,15, p=0.05, alternative="greater", correct=F)
@test netresult1[1] ≈ 0.8128703403598878   # p-value, from R
@test netresult1[3] == Dict("[0.05, 0.1)" => 2, "[0.1, 1.0]"  => 13)
@test netresult1[4][2] ≈ 68.03708830981597 # pseudo log-lik
@test netresult1[4][1] ≈ 29.34808731515701 atol=1e-5 # alpha
end

@testset "ticr! on data frame, network, minimum pvalue, dirichlet distribution" begin
net3 = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
# without optimizing branch lengths
netresult1 = ticr!(net3,d,false; quartetstat=:minpval);
@test netresult1[2] ≈ 8.589058727506838  # z stat
@test netresult1[1] ≈ 4.384234705965304e-18   # p-value
@test netresult1[3] == Dict("[0.05, 0.1)" => 1, "[0.0, 0.01)"  => 4, "[0.01, 0.05)" => 4, "[0.1, 1.0]"  => 6)
@test netresult1[4][2] ≈ 68.03708830981597 # pseudo log-lik
@test netresult1[4][1] ≈ 29.34808731515701 atol=1e-5 # alpha
end

end
