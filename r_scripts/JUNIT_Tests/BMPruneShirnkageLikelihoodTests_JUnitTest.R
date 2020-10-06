# This R script gives us the expected values for JUnit tests
# 'BMPruneShrinkageLikelihoodTest'
# (1) 'testBMPruneShrinkageLikelihood6Species5Traits'
# (2) 'testBMPruneShrinkageLikelihood3Species2Traits'
# (3) 'testBMPruneShrinkageLikelihood6Species5TraitsUltrametricTree'
# (4) 'testBMPruneShrinkageLikelihood3Species2TraitsPopSE'
# (5) 'testBMPruneWithWithoutShrinkageLikelihood'

library(mvMORPH)
library(PCMBase)
library(TreeSim)
library(mcmc3r)
library(corpcor)
library(phytools)

source('BMPruneLikelihoodJUnitTestUtils.R')

## 'BMPruneShrinkageLikelihoodTest'

main.path <- "./r_script/"
mcmcTree.path <- "/Applications/paml4.9j/bin"

### (1) testBMPruneShrinkageLikelihood6Species5Traits
# this tree only provides the topology and tip ages
# the internal node times will be extracted from the sampled tree in mcmcTree
tr <- read.tree(text = "((t4_2:1.467499928,t4_1:0.8):8.083170412,((t1_1:4.807085433,(t2_1:3.912659559,t3_1:3.912659559):0.8944258738):2.820091766,t6_1:3.0):1.923493141):1.705934219;")

trait.nr <- 5

dat <- matrix(c(-1.49423749961236, -1.01313598636176, 5.69171615806967, -1.76564797547353, -5.88260110808259,
             -3.25522309907749, -2.20443210635097, -1.45295055296761, -1.00502666401073, 2.73908099097237,
             -0.304175058530707, 5.68185395502633, 13.0073187571802, 8.20086495770802, 0.305337545316146,
             -10.2246819199159, -5.9742892461524, 1.38164443830931, -1.73557048290321, 0.459868047614318,
             0.330407167242628, 7.24123831044803, 0.674662257070058, -4.82948242309325, 0.362214907528917,
             0.494318850519285, 0.305460656456528, -0.143939688662232, -5.64650554205318, 4.23511851679411), nrow = 6, ncol = 5)

# optimal shrinkage parameter -> 0.963523553937718
delta <- estimate.lambda(x = dat, verbose = FALSE)

# trait correlations estimated by shrinlkage method
shrinkage.correlation.matrix <- diag(trait.nr) * delta + (1 - delta) * cor(dat)

# mcmcTree
file.name <- "testBMPruneShrinkageLikelihood6Species5Traits"

# creat an array for tip ages with species names which will be used to write input files of mcmcTree
# c("t4_2^9.55067034", "t4_1^8.883170412", "t1_1^9.55067034", "t2_1^9.55067034", "t3_1^9.55067034", "t6_1^4.923493141")
ages.str <- get.tip.ages.str(tr)

# prepare tree file, control file and alignment file for mcmcTree
write.script.for.mcmcTree(main.path, file.name, ages.str, NULL, dat, shrinkage.correlation.matrix, tr, "1 1", "10", "10", "1000", NULL)

# prepare bash script to run mcmcTree
write.bash.4.mcmcTree(mcmcTree.path, main.path, file.name)

# get results from the mcmc file
# node times (t_n) and branch rate (mu)
# Gen	t_n7	t_n8	t_n9	t_n10	t_n11	mu	lnL
# 10000	42.7712412	9.1379169	19.4908108	16.0044120	5.7809186	1.0479483	-76.655
# -> the expected likelihood is -76.655

# Note: there are 2 fossils, i.e. "t4_1 = 0.6674999280" and "t6_1 = 4.6271771990"
# check the tip ages by get.tip.ages(tr)
# which makes node times t4_1 = 9.1379169 - 0.6674999280 and t6_1 = 19.4908108 - 4.6271771990
# -> this is the tree to specify in the Java unit test
res.tree <- read.tree(text = "((t4_2:9.1379169,t4_1:8.470416972):33.6333243,((t1_1:16.0044120,(t2_1:5.7809186,t3_1:5.7809186):10.2234934):3.4863988,t6_1:14.863633601):23.2804304):0.0;")

# (2) 'testBMPruneShrinkageLikelihood3Species2Traits'
# this tree only provides the topology and tip ages
# the internal node times will be extracted from the sampled tree in mcmcTree
tr <- read.tree(text = "((A:8.362605913,B:8.362605913):5.94554475,C:2.510587814);")

trait.nr <- 2

dat <- matrix(c(1, 3, 2, 2, 5, 4), nrow = 3, ncol = 2)

# optimal shrinkage parameter -> 0.25925925925926
delta <- estimate.lambda(x = dat, verbose = FALSE)

# trait correlations estimated by shrinlkage method
shrinkage.correlation.matrix <- diag(trait.nr) * delta + (1 - delta) * cor(dat)

# mcmcTree
file.name <- "testBMPruneShrinkageLikelihood3Species2Traits"

# creat an array for tip ages with species names which will be used to write input files of mcmcTree
# c("A^14.308150663", "B^14.308150663", "C^2.510587814" )
ages.str <- get.tip.ages.str(tr)

# prepare tree file, control file and alignment file for mcmcTree
write.script.for.mcmcTree(main.path, file.name, ages.str, NULL, dat, shrinkage.correlation.matrix, tr, "1 1", "10", "10", "1000", NULL)

# prepare bash script to run mcmcTree
write.bash.4.mcmcTree(mcmcTree.path, main.path, file.name)

# get results from the mcmc file
# node times (t_n) and branch rate (mu)
# Gen	t_n4	t_n5	mu	lnL
# 10000	55.3678474	12.4420263	0.1543038	-8.128
# -> the expected likelihood is -8.128

# Note: C is a fossil
# so the node time should be (its parental node time - its tip date), i.e. (55.36784745 - 11.79756)
# check the tip age by get.tip.ages(tr)[3] -> 11.79756
# -> this is the tree to specify in the Java unit test
res.tree = read.tree(text = "((A:12.4420263,B:12.4420263):42.9258211,C:43.5702874);")

# (3) 'testBMPruneShrinkageLikelihood6Species5TraitsUltrametricTree'
# this tree only provides the topology and tip ages
# the internal node times will be extracted from the sampled tree in mcmcTree
tr <- read.tree(text = "((t4_2:1.467499928,t4_1:1.467499928):8.083170412,((t1_1:4.807085433,(t2_1:3.912659559,t3_1:3.912659559):0.8944258738):2.820091766,t6_1:7.627177199):1.923493141):0.0;")

trait.nr <- 5

dat <- matrix(c(-1.49423749961236, 1.61438393261534, 5.69171615806967, -1.76564797547353, -5.88260110808259, -11.0382115821499,
                -2.20443210635097, -0.600823508278549, -1.00502666401073, 2.73908099097237, -0.304175058530707, 7.74032190916081,
                13.0073187571802, 6.76701308252399, 0.305337545316146, -10.2246819199159, -5.9742892461524, 2.29305900150846,
                 -1.73557048290321, 1.35749484021919, 0.330407167242628, 7.24123831044803, 0.674662257070058, -7.36164535234136,
                 0.362214907528917, 0.0284817981254237, 0.305460656456528, -0.143939688662232, -5.64650554205318, 7.52792561895395), nrow = 6, ncol = 5)

# optimal shrinkage parameter -> 0.879396295242854
delta <- estimate.lambda(x = dat, verbose = FALSE)

# trait correlations estimated by shrinlkage method
shrinkage.correlation.matrix <- diag(trait.nr) * delta + (1 - delta) * cor(dat)

# mcmcTree
file.name <- "testBMPruneShrinkageLikelihood6Species5TraitsUltrametricTree"

# creat an array for tip ages with species names which will be used to write input files of mcmcTree
ages.list <- get.extant.str(tr)

# prepare tree file, control file and alignment file for mcmcTree
write.script.for.mcmcTree(main.path, file.name, NULL, ages.list, dat, shrinkage.correlation.matrix, tr, "0 0", "0", "1", "10", NULL)

# prepare bash script to run mcmcTree
write.bash.4.mcmcTree(mcmcTree.path, main.path, file.name)

# get results from the mcmc file
# node times (t_n) and branch rate (mu)
# Gen	t_n7	t_n8	t_n9	t_n10	t_n11	mu	lnL
# 10	55.7462536	4.1813368	24.0712897	16.4777971	2.4163095	0.0467689	-538.956
# -> the expected likelihood is -538.956
# -> this is the tree to specify in the Java unit test
res.tree = read.tree(text = "((t4_2:4.1813368,t4_1:4.1813368):51.5649168,((t1_1:16.4777971,(t2_1:2.4163095,t3_1:2.4163095):14.0614876):7.5934926,t6_1:24.0712897):31.6749639):0.0;")

# (4) 'testBMPruneShrinkageLikelihood3Species2TraitsPopSE'
tr <- read.tree(text = "((A:23.0058179,B:23.0058179):14.350951,C:37.3567689);")

trait.nr <- 2

dat <- matrix(c(-2.62762948691895, -1.50846427625826, -0.226074849617958,
                -1.56292164859448, -1.59482814741543, -2.11000367246907), nrow = 3, ncol = 2)

# optimal shrinkage parameter -> 0.303252482035699
delta <- estimate.lambda(x = dat, verbose = FALSE)

# trait correlations estimated by shrinlkage method
shrinkage.correlation.matrix <- diag(trait.nr) * delta + (1 - delta) * cor(dat)

# mcmcTree
file.name <- "testBMPruneShrinkageLikelihood3Species2TraitsPopSE"

# creat an array for tip ages with species names which will be used to write input files of mcmcTree
# c("A^37.3567689", "B^37.3567689", "C^37.3567689")
ages.str <- get.tip.ages.str(tr)
ages.list <- get.extant.str(tr)
# prepare tree file, control file and alignment file for mcmcTree
# population variance = 0.3
write.script.for.mcmcTree(main.path, file.name,  NULL, ages.list, dat, shrinkage.correlation.matrix, tr, "0 0", "0", "1", "10", 0.3)

# prepare bash script to run mcmcTree
write.bash.4.mcmcTree(mcmcTree.path, main.path, file.name)

# output from mcmcTree
# Gen	t_n4	      t_n5	      mu	         lnL
# 10	  39.3676799	23.3161955	0.3465491	-9.731
res.tree = read.tree(text = "((A:23.3161955,B:23.3161955):16.0514844,C:39.3676799);")
# Note that mu is the branch rate, which gives the following value in .java
# RealParameter colorValues = new RealParameter(new Double[] {0.3465491});

# (5) 'testBMPruneWithWithoutShrinkageLikelihood'
# this tree provides the topology and is used to simulate continuous traits
# the internal node times will be extracted from the sampled tree in mcmcTree
tr <- read.tree(text = "((t2:1.38278273,((t4:0.2275519885,t5:0.2275519885):0.2593961434,t3:0.4869481318):0.8958345985):0.307024997,t1:1.689807727);")

trait.nr <- 9

# define a BM model with 9 traits and 1 regime
modelBM <- PCM("BM", k = trait.nr, regimes = 1)
modelBM$Sigmae_x <- NULL
# specify the BM model parameters by randomly drawing from distributions 
set.seed(123)
modelBM$X0[1:trait.nr] = rnorm(trait.nr, mean = 0.0, sd = 2.0)
this.variance =  rlnorm(trait.nr, meanlog = 1.0, sdlog = 0.3)
this.covariance = rnorm(trait.nr * (trait.nr - 1) / 2, mean = 0.0, sd = 0.5)
# populate the Sigma matrix in PCMBase package
for (i in 1:trait.nr) {
  modelBM$Sigma_x[,,1][i,i] = this.variance[i]
}
k = 1
for (i in 1:(trait.nr-1)) {
  for (j in (i+1):trait.nr) {
    modelBM$Sigma_x[,,1][i,j] = this.covariance[k]
    k = k + 1
  }
}

# simulate traits
set.seed(123)
full.dat = PCMSim(tr, modelBM, X0 = modelBM$X0)
dat = t(full.dat[ , 1:PCMTreeNumTips(tr)])
# dat <- matrix(c(0.996206925293026, -2.65052276370749, -3.93588311090379, -5.59242339046684, -2.53243812537148,
#                 -2.91270268270845, -1.5792798388905, -6.4085861951559, -5.36037547409755, 5.72854673050833,
#                  1.35994187362101, 6.82985203851052, 9.63880711538996, 8.7939796440525, 5.53576146755777,
#                  0.102975767496277, 1.79880717599135, 1.31978890202522, 3.92494613908588, 1.16992957583461,
#                  -2.43052504782334, -0.867723116858107, -0.637507941646737, 0.0934316079593845, -0.550563880149175,
#                  2.98602562445202, 7.49844240046301, 8.33514902203189, 8.15176984162606, 2.86099545245125,
#                  -1.03975429308448, 4.50781790071978, 8.4164152934786, 6.70696946959889, 11.752265679971,
#                  -10.8901344038457, 0.520102801698764, 0.429643092693121, -2.11120755504519, -1.08417023425239,
#                  -1.30156469246422, -0.00971374138577277, -1.47520049052425, 0.74040432381528, -4.43344316554345), nrow = 5, ncol = 9)

# optimal shrinkage parameter -> 0.582249579199718
delta <- estimate.lambda(x = dat, verbose = FALSE)

# trait correlations estimated by shrinlkage method
shrinkage.correlation.matrix <- diag(trait.nr) * delta + (1 - delta) * cor(dat)

# mcmcTree
file.name <- "testBMPruneWithShrinkageLikelihood"

# creat an array for tip ages with species names which will be used to write input files of mcmcTree
ages.list <- get.extant.str(tr)

# prepare tree file, control file and alignment file for mcmcTree
write.script.for.mcmcTree(main.path, file.name, NULL, ages.list, dat, shrinkage.correlation.matrix, tr, "0 0", "0", "1", "10", NULL)

# prepare bash script to run mcmcTree
write.bash.4.mcmcTree(mcmcTree.path, main.path, file.name)

# get results from the mcmc file
#Gen	t_n6	t_n7	t_n8	t_n9	mu	lnL
#10	41.4589006	21.5238014	11.9789202	2.7950731	0.2125345	-101.547
# -> the expected likelihood is -101.547
res.tree1 = read.tree(text = "((t2:21.5238014,((t4:2.7950731,t5:2.7950731):9.1838471,t3:11.9789202):9.5448812):19.9350992,t1:41.4589006);")
# Note that mu is the branch rate, which gives the following value in .java
# RealParameter sigmasq1 = new RealParameter(new Double[] {0.2125345});

# if not using shrinkage, the trait rate matrix is the true values used when simulating the data 
true.correlation.matrix <- modelBM$Sigma_x[,,1] %*% t(modelBM$Sigma_x[,,1])

file.name <- "testBMPruneWithoutShrinkageLikelihood"

# prepare tree file, control file and alignment file for mcmcTree
write.script.for.mcmcTree(main.path, file.name, NULL, ages.list, dat, true.correlation.matrix, tr, "0 0", "0", "1", "10", NULL)

# prepare bash script to run mcmcTree
write.bash.4.mcmcTree(mcmcTree.path, main.path, file.name)

# output from mcmcTree
# Gen	t_n6	t_n7	t_n8	t_n9	mu	lnL
# 10	37.3524582	21.0771224	13.3011428	1.2011247	0.2061328	-102.655
res.tree2 = read.tree(text = "((t2:21.0771224,((t4:1.2011247,t5:1.2011247):12.1000181,t3:13.3011428):7.7759796):16.2753358,t1:37.3524582);")
# Note that mu is the branch rate, which gives the following value in .java
# RealParameter colorValues2 = new RealParameter(new Double[] {0.2061328});


