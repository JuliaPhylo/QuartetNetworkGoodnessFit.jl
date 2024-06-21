var documenterSearchIndex = {"docs":
[{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"using QuartetNetworkGoodnessFit\nENV[\"COLUMNS\"] = 85 # displaysize(stdout)[2] is 80 by default","category":"page"},{"location":"man/gof/#Goodness-of-fit-of-a-candidate-network-1","page":"goodness of fit","title":"Goodness of fit of a candidate network","text":"","category":"section"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"The example data set in this package is very small (5 taxa) on purpose: to make this tutorial run quickly. The data list all four-taxon sets, and the concordance factor (CF) for each of the 3 topologies on these 4 taxa, that is, the proportion of genes estimated to have each 4-taxon unrooted topology.","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"using QuartetNetworkGoodnessFit, CSV\nqCF = CSV.read(joinpath(dirname(pathof(QuartetNetworkGoodnessFit)), \"..\",\"test\",\"example_qCF_5taxa.csv\"));\nqCF","category":"page"},{"location":"man/gof/#testing-candidate-networks-1","page":"goodness of fit","title":"testing candidate networks","text":"","category":"section"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"Let's say we have two candidate networks to test: one that is just a tree (net0), and the other that has one reticulation (net1).","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"using PhyloNetworks\nnet1 = readTopology(\"((((D,C),((A,B))#H1),(#H1,E)),O);\");\nnet0 = readTopology(\"((((D,C),(A,B)),E),O);\"); # also: majorTree(net1)","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"net0 has clades AB and CD sister to each other, with E sister to ABCD and an outgroup O. net1 has net0 as its major tree, and has a reticulation going from E into the stem lineage to AB.","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"The goal is to test the goodness-of-fit of each topology.","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"It's best to start with a small number of simulations for the test, to check that there are no issue. Then we can re-run with a larger number. Note that the first use of a function includes time to compile the function, so it takes longer than the second use. Below, the option true means that we want to optimize the branch lengths in the network, in coalescent units, before quantifying the goodness-of-fit.","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"res0 = quarnetGoFtest!(net0, qCF, true; seed=234, nsim=2);\nnothing # hide","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"Before re-running with more simulations, we can update the network net0 to the version that has optimized branch lengths:","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"net0 = res0[5]","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"Now we re-run the test using the option false to not re-optimize branch lengths. We use nsim=200 simulations below to make this example faster. For a real data analysis, delete the nsim option to use the default instead (1000) or specify higher value.","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"res0 = quarnetGoFtest!(net0, qCF, false; seed=234, nsim=200);","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"In the result, the first 3 numbers are the p-value of the overall goodness-of-fit test, the uncorrected z-value, and the estimate of σ to correct the z-value for dependence:","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"res0[[1,2,3]] # p-value, uncorrected z, σ","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"Here, we have evidence that the tree net0 does not fit the data adequately. Note that σ=1 corresponds to independent quartet outlier p-values. Typical estimates of σ are quite larger than 1, increase with the number of taxa, and increase with longer branch lengths (for a given number of taxa).","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"The next element is the list of outlier p-values, which was also added to the data frame that contains our data:","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"res0[4]\nqCF[:,[:t1,:t2,:t3,:t4,:p_value]]","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"The very small overall p-value indicated an excess of outlier four-taxon sets. Here we see what these outliers are: all four-taxon sets containing O, E and either A or B are strong outliers (very small outlier p-values).","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"In fact, the CF data were simulated on net1, in which the AB clade received gene flow from E. (See below if interested in the simulation). We can re-run an analysis using net1 this time:","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"res1 = quarnetGoFtest!(net1, qCF, true; seed=721, nsim=200);\nres1[[1,2,3]]","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"This network is found to provide an adequate absolute fit (as expected from how the data were simulated).","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"Note that after optimization of branch lengths and γs to best fit the CF data, the network was \"ultrametrized\" along its major tree to adhere to the hybrid-lambda simulator's requirement that the network is ultrametric:","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"res1[5]","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"For more options, see quarnetGoFtest!, such as for the outlier test statistic (G or likelihood ratio test by default).","category":"page"},{"location":"man/gof/#parallel-computations-1","page":"goodness of fit","title":"parallel computations","text":"","category":"section"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"For larger networks, the large number of four-taxon sets causes the test to run more slowly. The computations can be parallelized by providing more processors to Julia. For instance, this can be done by starting julia with the -p option:","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"julia -p 3 # 3 worker processors, 4 processors total","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"To check for progress during the test, we can check for a new directory with a name starting with jl_, such as jl_0CVOfE (the end of this name is randomly generated). Then we can check the files in this directory: their names are indicative of the simulation replicate number that is currently being processed. In the example below, there are 3 files in my jl_0CVOfE directory (because of using 3 workers), and these files show processing replicate numbers 106, 440 and 773 (out of the 1000 replicates by default). Each contains 200 lines (because it contains 1 gene tree per line, since the original data had 200 genes). It means that the simulation-based test is at about 1/3 of its way. By the way, these files should not be modified. At the end of the simulation-based test, these temporary files and folder are automatically deleted.","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"$ wc -l jl_0CVOfE/*\n     200 jl_0CVOfE/genetrees_rep106_coal_unit\n     200 jl_0CVOfE/genetrees_rep440_coal_unit\n     200 jl_0CVOfE/genetrees_rep773_coal_unit","category":"page"},{"location":"man/gof/#quartet-concordance-factor-simulation-1","page":"goodness of fit","title":"quartet concordance factor simulation","text":"","category":"section"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"The data in qCF were simulated using hybrid-Lambda, which comes with QuartetNetworkGoodnessFit. In the code below, hybrid-lambda is called within julia, asking for the simulation of 200 genes (-num 200).","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"Note that hybrid-lambda assumes that the tree or network is ultrametric, but runs even if this assumption is not met.","category":"page"},{"location":"man/gof/#","page":"goodness of fit","title":"goodness of fit","text":"net1 = readTopology(\"((((D:0.2,C:0.2):2.4,((A:0.4,B:0.4):1.1)#H1:1.1::0.7):2.0,(#H1:0.0::0.3,E:1.5):3.1):1.0,O:5.6);\");\nhl = QuartetNetworkGoodnessFit.hybridlambda # path to hybrid-lambda simulator, on local machine\nnet1HL = hybridlambdaformat(net1) # format for the network, to use below by hybrid-lambda\nrun(`$hl -spcu \"((((D:0.2,C:0.2)I1:2.4,((A:0.4,B:0.4)I2:1.1)H1#0.7:1.1)I3:2.0,(H1#0.7:0.0,E:1.5)I4:3.1)I5:1.0,O:5.6)I6;\" -num 200 -seed 123 -o \"genetrees\"`)\ntreelist = readMultiTopology(\"genetrees_coal_unit\")\nobsCF = writeTableCF(countquartetsintrees(treelist)...)","category":"page"},{"location":"lib/public/#public-documentation-1","page":"public","title":"public documentation","text":"","category":"section"},{"location":"lib/public/#","page":"public","title":"public","text":"Documentation for QuartetNetworkGoodnessFit's public (exported) functions.","category":"page"},{"location":"lib/public/#index-1","page":"public","title":"index","text":"","category":"section"},{"location":"lib/public/#","page":"public","title":"public","text":"Pages = [\"public.md\"]","category":"page"},{"location":"lib/public/#functions-1","page":"public","title":"functions","text":"","category":"section"},{"location":"lib/public/#","page":"public","title":"public","text":"Modules = [QuartetNetworkGoodnessFit]\nPrivate = false\nOrder   = [:function]","category":"page"},{"location":"lib/public/#QuartetNetworkGoodnessFit.quarnetGoFtest!-Tuple{HybridNetwork,DataFrames.DataFrame,Bool}","page":"public","title":"QuartetNetworkGoodnessFit.quarnetGoFtest!","text":"quarnetGoFtest!(net::HybridNetwork, df::DataFrame, optbl::Bool; quartetstat=:LRT, correction=:simulation, seed=1234, nsim=1000, verbose=false, keepfiles=false)\nquarnetGoFtest!(net::HybridNetwork, dcf::DataCF,   optbl::Bool; kwargs...)\n\nGoodness-of-fit test for the adequacy of the multispecies network coalescent, to see if a given population or species network explains the quartet concordance factor data adequately. The network needs to be of level 1 at most (trees fullfil this condition), and have branch lengths in coalescent units. The test assumes a multinomial distribution for the observed quartet concordance factors (CF), such that information on the number of genes for each four-taxon set (ngenes field) must be present.\n\nFor each four-taxon set, an outlier p-value is obtained by comparing a test statistic (-2log likelihood ratio by default) to a chi-square distribution with 2 degrees of freedom (see below for other options).\n\nThe four-taxon sets are then categorized as outliers or not, according to their outlier p-values (outlier if p<0.05). Finally, a one-sided goodness-of-fit test is performed on the frequency of outlier 4-taxon sets to see if there are more outliers than expected. The z-value for this test corresponds to the null hypothesis that 5% outlier p-values are < 0.05 (versus more than 5%):\n\nz = fracmathrmproportionoutliers - 005sqrt005 times 095mathrmnumber4taxonsets\n\nThis z-value corresponds to a test that assumes independent outlier p-values across 4-taxon sets: it makes no correction for dependence.\n\nTo correct for dependence with correction=:simulation, the distribution of z-values is obtained by simulating gene trees under the coalescent along the network (after branch length optimization if optbl=true) using hybrid-Lambda. The z-score is calculated on each simulated data set. Under independence, these z-scores have mean 0 and variance 1. Under dependence, these z-scores still have mean 0, but an inflated variance. This variance σ² is estimated from the simulations, and the corrected p-value is obtained by comparing the original z value to N(0,σ²). When correction=:none, σ is taken to be 1 (independence): not recommended!\n\nThe first version takes a DataFrame where each row corresponds to a given four-taxon set. The data frame is modified by having an additional another column containing the p-values corresponding to each four-taxon set.\nThe second version takes a DataCF object and modifies it by updating the expected concordance factors stored in that object.\n\nNote that net is not modified.\n\narguments\n\noptbl: when false, branch lengths in net are taken as is, and need to be in coalescent units. When optbl=true, branch lengths in net are optimized, to optimize the pseudo log likelihood score as in SNaQ (see here). In both cases, any missing branch length is assigned a value with ultrametrize!, to make the major tree ultrametric, in an attempt to adhere to the current requirements of the hybrid-lambda simulator. Missing branch lengths may arise if they are not identifiable, such as lengths of external branches if there is a single allele per taxon. The network is returned as part of the output.\n\nkeyword arguments\n\nquartetstat: the test statistic used to obtain an outlier p-value for each four-taxon set, which is then compared to a chi-squared distribution with 2 degrees of freedom to get a p-value. The default is :LRT for the likelihood ratio: 2n_mathrmgenes sum_j=1^3 hat p_j (loghat p_j - log p_j) where p_j is the quartet CF expected from the network, and hat p_j is the quartet CF observed in the data.   Alternatives are :Qlog for the Qlog statistics (Lorenzen, 1995): 2n_mathrmgenes sum_j=1^3 frac(hat p_j - p_j)^2p_j (loghat p_j - log p_j)   and :pearson for Pearon's chi-squared statistic, which behaves poorly when the expected count is low (e.g. less than 5): n_mathrmgenes sum_j=1^3 frac(hat p_j - p_j)^2 p_j\ncorrection=:simulation to correct for dependence across 4-taxon. Use :none to turn off simulations and the correction for dependence.\nseed=1234: master seed to control the seeds for gene tree simulations.\nnsim=1000: number of simulated data sets. Each data set is simulated to have the median number of genes that each 4-taxon sets has data for.\nverbose=false: output from hybrid-Lambda and other output is suppressed. Turn to true to diagnose potential issues.\nkeepfiles=false: files generated by hybrid-Lambda are stored in a temporary folder, whose name starts with jl_ and placed in the current directory. Each of the 1000 simulations creates a file, that is deleted when is no longer necessary. Turn keepfiles=true to keep these gene tree files.\n\noutput\n\np-value of the overall goodness-of-fit test (corrected for dependence if requested)\nuncorrected z value test statistic\nestimated σ for the test statistic used for the correction (1.0 if no correction)\na vector of outlier p-values, one for each four-taxon set\nnetwork (first and second versions): net with loglik field updated if optbl is false;  copy of net with optimized branch lengths and loglik if optbl is true\nin case correction = :simulation, vector of simulated z values (nothing if correction = :none). These z-values could be used to calculate an empirical p-value (instead of the p-value in #1), as the proportion of simulated z-values that are ⩾ the observed z-value (in #2).\n\nreferences\n\nRuoyi Cai & Cécile Ané (in prep). Assessing the fit of the multi-species network coalescent to multi-locus data.\nLorenzen (1995). A new family of goodness-of-fit statistics for discrete multivariate data. Statistics & Probability Letters, 25(4):301-307. doi: 10.1016/0167-7152(94)00234-8\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#QuartetNetworkGoodnessFit.ticr!-Tuple{HybridNetwork,DataFrames.DataFrame,Bool}","page":"public","title":"QuartetNetworkGoodnessFit.ticr!","text":"ticr!(net::HybridNetwork, df::DataFrame, optbl::Bool; quartetstat, test)\nticr!(net::HybridNetwork, dcf::DataCF,   optbl::Bool; quartetstat, test)\nticr(dcf::DataCF, quartetstat::Symbol, test::Symbol)\n\nGoodness-of-fit test for the adequacy of the multispecies network coalescent, to see if a given population or species network explains the quartet concordance factor data adequately. See Stenz et al 2015 and addendum for the method on trees, from which the acronym TICR was derived: Tree Incongruence Checking with R.\n\nThe tree / network needs to have branch lengths in coalescent units, and must be of level 1. It must be fully resolved, such that the \"major\" quartet is clearly defined, even though branch lengths can be 0.\n\nThe Dirichlet distribution is used to fit the distribution of the observed quartet concordance factors (CF), with a concentration parameter estimated from the data.\n\nThe four-taxon sets are then binned into categories according to their outlier p-values: 0-0.01, 0.01-0.05, 0.05-0.10, and 0.10-1. Finally, a one-sided goodness-of-fit test is performed on these binned frequencies to determine if they depart from an expected 5% p-values being below 0.05 (versus more than 5%), using the default test = :onesided. Alternatively, the option test = :goodness carries out the overall goodness-of-fit chi-square test proposed by Stenz et al. (2015), to test for any kind of departure from the expected frequencies (0.01, 0.04, 0.05, 0.90) across all 4 bins.\n\nThe first version takes a DataFrame where each row corresponds to a given four-taxon set. The data frame is modified by having an additional another column containing the outlier p-values corresponding to each four-taxon set.\nThe second version takes a DataCF object and modifies it by updating the expected concordance factors stored in that object.\nThe last version (which all others call) assumes that the expected concordance factors in the DataCF object are correctly calculated from the test network.\n\narguments\n\noptbl: when false, the loglik field of net is updated but branch lengths are taken as is. When true, a copy of net is used to conduce the test, with updated branch lengths (in coalescent units) and updated loglik. This network is returned.\nquartetstat = :maxCF: test statistic used to obtain an outlier p-value for each four-taxon set. By default, it is the absolute difference between the observed CF and expected CF of the major resolution (which has the largest CF) if quartetstat=:maxCF, as used in Stenz et al. (2015). The other option is quartetstat=:minpval, in which case the absolute difference between the observed CF and expected CF is calculated for all 3 resolutions, leading to 3 (non-independent) p-values. The outlier p-value for a given four-taxon set is taken as the minimum of these 3 p-values. This option can detect all kinds of departure from the network model for a given four-taxon set, but is not recommended because it is very liberal.\ntest = :onesided: the overall goodness-of-fit test performed on binned frequencies of quartet outlier p-values; see above.\n\noutput\n\np-value of the overall goodness-of-fit test\ntest statistic (z value)\na dictionary of the count of p-values in each of the four category\na vector of two test parameters for the Dirichlet distribution:\nvalue of the concentration parameter α\nvalue of the pseudo likelihood (optimized at α)\na vector of outlier p-values, one for each four-taxon set\nnetwork (first and second versions): net with loglik field updated if optbl is false;  copy of net with optimized branch lengths and loglik if optbl is true\n\nreferences\n\nNWM Stenz, B Larget, DA Baum and C Ané (2015). Exploring tree-like and non-tree-like patterns using genome sequences: An example using the inbreeding plant species Arabidopsis thaliana (L.) Heynh. Systematic Biology, 64(5):809-823. doi: 10.1093/sysbio/syv039\n\nsee also this addendum\n\n\n\n\n\n","category":"method"},{"location":"#QuartetNetworkGoodnessFit.jl-1","page":"home","title":"QuartetNetworkGoodnessFit.jl","text":"","category":"section"},{"location":"#","page":"home","title":"home","text":"Julia package to measure the goodness of fit of a phylogenetic network to data on subsets of 4 tips. It depends on the PhyloNetworks package.","category":"page"},{"location":"#","page":"home","title":"home","text":"For a tutorial, see the manual:","category":"page"},{"location":"#","page":"home","title":"home","text":"Pages = [\n    \"man/gof.md\",\n]\nDepth = 1","category":"page"},{"location":"#","page":"home","title":"home","text":"References:","category":"page"},{"location":"#","page":"home","title":"home","text":"Ruoyi Cai & Cécile Ané (in prep). Assessing the fit of the multi-species network coalescent to multi-locus data.\nNoah W. M. Stenz, Bret Larget, David A. Baum, and Cécile Ané (2015). Exploring tree-like and non-tree-like patterns using genome sequences: An example using the inbreeding plant species Arabidopsis thaliana (L.) Heynh. Systematic Biology, 64(5):809-823.\nAddendum describing a modification to the model in the original TICR test.","category":"page"},{"location":"#functions-1","page":"home","title":"functions","text":"","category":"section"},{"location":"#","page":"home","title":"home","text":"Pages = [\"lib/public.md\", \"lib/internals.md\"]\nOrder = [:function]","category":"page"},{"location":"lib/internal/#internal-documentation-1","page":"internal","title":"internal documentation","text":"","category":"section"},{"location":"lib/internal/#","page":"internal","title":"internal","text":"Documentation for QuartetNetworkGoodnessFit's internal functions. Those functions are not exported, but can still be used (like: QuartetNetworkGoodnessFit.foo() for a function named foo()).","category":"page"},{"location":"lib/internal/#index-1","page":"internal","title":"index","text":"","category":"section"},{"location":"lib/internal/#","page":"internal","title":"internal","text":"Pages = [\"internals.md\"]","category":"page"},{"location":"lib/internal/#functions-1","page":"internal","title":"functions","text":"","category":"section"},{"location":"lib/internal/#","page":"internal","title":"internal","text":"Modules = [QuartetNetworkGoodnessFit]\nPublic  = false\nOrder   = [:function]","category":"page"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.dirichlet_max-Tuple{DataCF}","page":"internal","title":"QuartetNetworkGoodnessFit.dirichlet_max","text":"dirichlet_max(dcf::DataCF)\n\nCalculate outlier p-values, one for each four-taxon set, using the maximum concordance factor under a Dirichlet distribution. Used by ticr!.\n\noutput\n\nvector of outlier p-values, one for each 4-taxon set\nvalue of the concentration parameter α\nvalue of the pseudo likelihood (optimized at α)\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.dirichlet_min-Tuple{DataCF}","page":"internal","title":"QuartetNetworkGoodnessFit.dirichlet_min","text":"dirichlet_min(dcf::DataCF)\n\nFirst calculate outlier p-values using each of the three concordance factors for each quartet under a Dirichlet distribution, then take the smallest p-value among the three, as the outlier p-value for each four-taxon set.\n\noutput\n\na vector of outlier p-values, one for each quartet\nvalue of the concentration parameter α\nvalue of the pseudo likelihood (optimized at α)\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.expectedCF_ordered","page":"internal","title":"QuartetNetworkGoodnessFit.expectedCF_ordered","text":"expectedCF_ordered(dcf::DataCF, net::HybridNetwork, suffix=\"\"::AbstractString)\n\nExpected quartet concordance factors in dcf, but ordered as they would be if output by PhyloNetworks.countquartetsintrees. Output:\n\n2-dimentional SharedArray (number of 4-taxon sets x 3). dcf.quartet[i].qnet.expCF[j] for 4-taxon set i and resolution j is stored in row qi and column k if qi is the rank of 4-taxon set i (see PhyloNetworks.quartetrank). This rank depends on how taxa are ordered.\nvector of taxon names, whose order matters. These are tip labels in net with suffix suffix added, then ordered alphabetically, or numerically if taxon names can be parsed as integers.\n\n\n\n\n\n","category":"function"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.multinom_lrt!-Tuple{AbstractArray{Float64,1},Array{Quartet,1}}","page":"internal","title":"QuartetNetworkGoodnessFit.multinom_lrt!","text":"multinom_lrt!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})\nmultinom_lrt!(pval::AbstractVector{Float64}, obsCF, expCF::AbstractMatrix{Float64})\n\nCalculate outlier p-values (one per four-taxon set) using the likelihood ratio test under a multinomial distribution for the observed concordance factors.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.multinom_pearson!-Tuple{AbstractArray{Float64,1},Array{Quartet,1}}","page":"internal","title":"QuartetNetworkGoodnessFit.multinom_pearson!","text":"multinom_pearson!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})\nmultinom_pearson!(pval::AbstractVector{Float64}, obsCF, expCF::AbstractMatrix{Float64})\n\nCalculate outlier p-values (one per four-taxon set) using Pearson's chi-squared statistic under a multinomial distribution for the observed concordance factors.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.multinom_qlog!-Tuple{AbstractArray{Float64,1},Array{Quartet,1}}","page":"internal","title":"QuartetNetworkGoodnessFit.multinom_qlog!","text":"multinom_qlog!(pval::AbstractVector{Float64}, quartet::Vector{Quartet})\nmultinom_qlog!(pval::AbstractVector{Float64}, obsCF, expCF::AbstractMatrix{Float64})\n\nCalculate outlier p-values (one per four-taxon set) using the Qlog statistic (Lorenzen, 1995), under a multinomial distribution for the observed concordance factors.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.quarnetGoFtest-Tuple{Array{Quartet,1},Function}","page":"internal","title":"QuartetNetworkGoodnessFit.quarnetGoFtest","text":"quarnetGoFtest(quartet::Vector{Quartet}, outlierp_fun!::Function)\nquarnetGoFtest(outlier_pvalues::AbstractVector)\n\nCalculate an outlier p-value for each quartet according to function outlierp_fun! (or take outlier-values as input: second version) and calculate the z-value to test the null hypothesis that 5% of the p-values are < 0.05, versus the one-sided alternative of more outliers than expected.\n\nSee quarnetGoFtest! for more details.\n\nOutput:\n\nz-value\noutlier p-values (first version only)\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.quarnetGoFtest_simulation-Tuple{HybridNetwork,DataCF,Function,Int64,Int64,Bool,Bool}","page":"internal","title":"QuartetNetworkGoodnessFit.quarnetGoFtest_simulation","text":"quarnetGoFtest_simulation(net::HybridNetwork, dcf::DataCF, outlierp_fun!::Function,\n                          seed::Int, nsim::Int, verbose::Bool, keepfiles::Bool)\n\nRun external program hybrid-Lambda (downloaded at installation, with path stored in QuartetNetworkGoodnessFit.hybridlambda) to simulate gene trees under the multispecies coalescent model along network net. The quartet concordance factors (CFs) from these simulated gene trees are used as input to outlierp_fun! to categorize each 4-taxon set as an outlier (p-value < 0.05) or not. For each simulated data set, a goodness-of-fit z-value is calculated by comparing the proportion of outlier 4-taxon sets to 0.05. The standard deviation of these z-values (assuming a mean of 0), and the z-values themselves are returned.\n\nUsed by quarnetGoFtest!.\n\nWarning: The quartet CFs expected from net are assumed to be stored in dcf.quartet[i].qnet.expCF. This is not checked.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.reroot!-Tuple{Any,Any}","page":"internal","title":"QuartetNetworkGoodnessFit.reroot!","text":"reroot!(net, refnet)\n\nReroot net to minimize the hardwired cluster distance between the net (with the new root position) and the reference network refnet. Candidate root positions are limited to internal nodes (excluding leaves) that are compatible with the direction of hybrid edges.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.ticr_optimalpha-Tuple{DataCF}","page":"internal","title":"QuartetNetworkGoodnessFit.ticr_optimalpha","text":"ticr_optimalpha(dcf::DataCF)\n\nFind the concentration parameter α by maximizing the pseudo-log-likelihood of observed quartet concordance factors. The model assumes a Dirichlet distribution with mean equal to the expected concordance factors calculated from a phylogenetic network (under ILS). These expected CFs are assumed to be already calculated, and stored in dcf.\n\nWhen calculating the pseudo-log-likelihood, this function checks the observed concordance factors for any values equal to zero: they cause a problem because the Dirichlet density is 0 at 0 (for concentration α > 1). Those 0.0 observed CF values are re-set to the minimum of:\n\nthe minimum of all expected concordance factors, and\nthe minimum of all nonzero observed concordance factors.\n\noutput\n\nmaximized pseudo-loglikelihood\nvalue of α where the pseudo-loglikelihood is maximized\nreturn code of the optimization\n\nThe optimization uses NLOpt, with the :LN_BOBYQA method. Optional arguments can tune the optimization differently: nloptmethod, xtol_rel (1e-6 by default), starting α value x_start (1.0 by default).\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/#QuartetNetworkGoodnessFit.ultrametrize!-Tuple{HybridNetwork,Bool}","page":"internal","title":"QuartetNetworkGoodnessFit.ultrametrize!","text":"ultrametrize!(net::HybridNetwork, verbose::Bool)\n\nAssign values to missing branch lengths in net to make the network time-consistent (all paths from the root to a given hybrid node have the same length) and ultrametric (all paths from the root to the tips have the same length), if possible.\n\nWarnings are given if it's not possible and if verbose is true.\n\nThe major tree is used to assign node ages. If an edge has a missing length, this length is changed to the following:\n\n0 if the edge is internal and non-minor,\nthe value needed to make the network time-consistent if the edge is minor,\nthe smallest value possible to make the network ultrametric if the edge is external.\n\n\n\n\n\n","category":"method"}]
}
