# author: Pau Bravo
# This R script gives us the expected values for JUnit tests
# 'OUUtilsComputeOUTMatOneTraitTest'
# (1) 'testComputeOUTMatSmallTree'
# (2) 'testComputeOUTMatSmallTreeNonUltra'
# (3) 'testComputeOUTMatLargeTreeNonUltra'

library(mvMORPH)

## The following three commands are aimed to modify mvMORPH package in order for it to print
## the final covariance matrix (multiplied by sigma^2) and weight matrix of the OU model to compare our implementation with mvMORPH's.
envirMORPH <- environment(mvMORPH::mvOU)
source(file="../mvOU.R")
environment(mvOU) <- envirMORPH

# The following command loads prior functions to calculate the covariance matrix,
# the weight matrix and the likelihood of the Hansen model to compare their results with mvMORPH's results
source(file = "../OUfunctions.R")

EPSILON <- 1e-04 # Tolerance value for comparing mvMORPH with PAU's functions

## OUMVNLikelihoodOneTraitFromUtilsTest

# ----- START: MVNUtils.computeMVNLk validation ----- #
# JUnit: OUOneTraitcomputeMVNLkTest

## R: rootIsRandVar=true (vcv="randomRoot")
## F: rootIsRandVar=false (vcv="fixedRoot")
## I: useRootMetaData=true (root=T)
## M: useRootMetaData=false (root=F)

### (1) testComputeOULikSmallTreeThreeOpt

tr <- paintSubTree(read.tree(text="(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);"), node=6, state="group1", anc.state="group0")
tr <- paintSubTree(tr, node=3, state="group2", anc.state="group0", stem=T)
dat <- c(4.1, 4.5, 5.9, 0.0)
plotSimmap(tr)

fitFI <- mvOU(tr, dat, model="OUM", param=list(root=T, vcv="fixedRoot"))
fitFM <- mvOU(tr, dat, model="OUM", param=list(root=F, vcv="fixedRoot"))
fitRI <- mvOU(tr, dat, model="OUM", param=list(root=T, vcv="randomRoot"))
fitRM <- mvOU(tr, dat, model="OUM", param=list(root=F, vcv="randomRoot"))

#### Values to use on JUNIT Tests
fitFI$alpha   # 31.20814
fitFI$sigma   # 1.248328
fitFI$theta   # c(theta0, group0, group1, group2) = c(2.228585e-40, -4.047373e-16, 4.300000e+00, 5.900000e+00)

fitFM$alpha   # 31.20814
fitFM$sigma   #  1.248328
fitFM$theta   # c(group0, group1, group2) = c(8.128044e-27, 4.300000e+00, 5.900000e+00)

fitRI$alpha   # 43.35287
fitRI$sigma   # 1.734117
fitRI$theta   # c(theta0, group0, group1, group2) = c(3.348633e-56, -1.903330e-16, 4.300000e+00, 5.900000e+00)

fitRM$alpha   # 43.35287
fitRM$sigma   # 1.734117
fitRM$theta   # c(group0, group1, group2) = c(2.297268e-37, 4.300000e+00, 5.900000e+00)

#### Values to compare with Java classe# Pau's functions output
##### Pau's functions output
regimes <- c(1, 1, 2, 0, 0, 0, 1); regimes <- factor(regimes);

ouCovFI <- varOU(ultTree, as.numeric(fitFI$alpha), T)
ouWFI <- weightMat(ultTree, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)

ouCovFM <- varOU(ultTree, as.numeric(fitFM$alpha), T)
ouWFM <- weightMat(ultTree, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)

ouCovRI <- varOU(ultTree, as.numeric(fitRI$alpha), F)
ouWRI <- weightMat(ultTree, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)

ouCovRM <- varOU(ultTree, as.numeric(fitRM$alpha), F)
ouWRM <- weightMat(ultTree, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

likFI <- computeLk(ult.data, ouWFI %*% fitFI$theta, ouCovFI, as.numeric(fitFI$sigma))
likFM <- computeLk(ult.data, ouWFM %*% fitFM$theta, ouCovFM, as.numeric(fitFM$sigma))
likRI <- computeLk(ult.data, ouWRI %*% fitRI$theta, ouCovRI, as.numeric(fitRI$sigma))
likRM <- computeLk(ult.data, ouWRM %*% fitRM$theta, ouCovRM, as.numeric(fitRM$sigma))

##### LogLikelihoods
fitFI$LogLik  # 2.148292
fitFM$LogLik  # 2.148292
fitRI$LogLik  # 2.148292
fitRM$LogLik  # 2.148292

##### Comparing mvMORPH output with Pau's R functions
all.equal(as.numeric(log(likFI)), as.numeric(fitFI$LogLik), tolerance=EPSILON)
all.equal(as.numeric(log(likFM)), as.numeric(fitFM$LogLik), tolerance=EPSILON)
all.equal(as.numeric(log(likRI)), as.numeric(fitRI$LogLik), tolerance=EPSILON)
all.equal(as.numeric(log(likRM)), as.numeric(fitRM$LogLik), tolerance=EPSILON)

### (2) testComputeOULikSmallTreeNonUltraOneOpt

tr2 <- read.tree(text="(((sp1:2.0, sp2:1.0):1.0, sp3:4.0):1.0, sp4:3.0);")
dat <- c(4.1, 4.5, 5.9, 0.0)
plot(tr2)

fitFI <- mvOU(tr2, dat, model="OU1", param=list(root=T, vcv="fixedRoot"))
fitFM <- mvOU(tr2, dat, model="OU1", param=list(root=F, vcv="fixedRoot"))
fitRI <- mvOU(tr2, dat, model="OU1", param=list(root=T, vcv="randomRoot"))
fitRM <- mvOU(tr2, dat, model="OU1", param=list(root=F, vcv="randomRoot"))

#### Values to use on JUNIT Tests
fitFI$alpha   # 0.7465763
fitFI$sigma   # 4.003551
fitFI$theta   # c(theta0, group0) = c(-33.591241, 6.449917)

fitFM$alpha   # 1.40338e-08
fitFM$sigma   # 1.237864
fitFM$theta   # 2.792045

fitRI$alpha   # 0.8609833
fitRI$sigma   # 4.601164
fitRI$theta   # c(theta0, group0) = c(-46.965464, 6.201598)

fitRM$alpha   # 0.7085376
fitRM$sigma   # 6.841867
fitRM$theta   # 3.586504

#### Values to compare with Java classe# Pau's functions output
##### Pau's functions output
regimes <- c(0, 0, 0, 0, 0, 0, 0); regimes <- factor(regimes);

ouCovFI <- varOU(tr2, as.numeric(fitFI$alpha), T)
ouWFI <- weightMat(tr2, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)

ouCovFM <- varOU(tr2, as.numeric(fitFM$alpha), T)
ouWFM <- weightMat(tr2, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)

ouCovRI <- varOU(tr2, as.numeric(fitRI$alpha), F)
ouWRI <- weightMat(tr2, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)

ouCovRM <- varOU(tr2, as.numeric(fitRM$alpha), F)
ouWRM <- weightMat(tr2, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

likFI <- computeLk(dat, ouWFI %*% fitFI$theta, ouCovFI, as.numeric(fitFI$sigma))
likFM <- computeLk(dat, ouWFM %*% fitFM$theta, ouCovFM, as.numeric(fitFM$sigma))
likRI <- computeLk(dat, ouWRI %*% fitRI$theta, ouCovRI, as.numeric(fitRI$sigma))
likRM <- computeLk(dat, ouWRM %*% fitRM$theta, ouCovRM, as.numeric(fitRM$sigma))

##### LogLikelihoods
fitFI$LogLik  # -7.630117
fitFM$LogLik  # -8.457486
fitRI$LogLik  # -7.63854
fitRM$LogLik  # -8.817273

##### Comparing mvMORPH output with Pau's R functions
all.equal(as.numeric(log(likFI)), as.numeric(fitFI$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likFM)), as.numeric(fitFM$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likRI)), as.numeric(fitRI$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likRM)), as.numeric(fitRM$LogLik), tolerance = EPSILON)

### (3) testComputeOULikLargeTreeNonUltraFiveOpt

tr3 <- paintSubTree(read.tree(text="(((((sp1:1.0,sp2:1.0):1.0,sp3:1.0):2.0,(sp4:1.0,sp5:1.0):3.0):2.0,sp6:6.0):1.0,sp7:3.0);"), node=12, state="group1", anc.state="group0", stem=T)
tr3 <- paintSubTree(tr3, node=12, state="group1", anc.state="group0", stem=T)
tr3 <- paintSubTree(tr3, node=13, state="group2", anc.state="group0", stem=T)
tr3 <- paintSubTree(tr3, node=6, state="group3", anc.state="group0", stem=T)
tr3 <- paintSubTree(tr3, node=7, state="group4", anc.state="group0", stem=T)
dat <- c(4.1, 4.5, 5.9, 0.0, 3.2, 2.5, 5.4)
plotSimmap(tr3)

fitFI <- mvOU(tr3, dat, model="OUM", param=list(root=T, vcv="fixedRoot"))
fitFM <- mvOU(tr3, dat, model="OUM", param=list(root=F, vcv="fixedRoot"))
fitRI <- mvOU(tr3, dat, model="OUM", param=list(root=T, vcv="randomRoot"))
fitRM <- mvOU(tr3, dat, model="OUM", param=list(root=F, vcv="randomRoot"))

fitFI$alpha   # 7.142986
fitFI$sigma   # 10.61289
fitFI$theta   # c(2.666338e-09, 5.900000e+00, 4.299999e+00, 1.600000e+00, 2.500000e+00, 5.400000e+00)

fitFM$alpha   # 7.142986
fitFM$sigma   # 10.61289
fitFM$theta   # c(5.900000e+00, 4.299999e+00, 1.600000e+00, 2.500000e+00, 5.400000e+00)

fitRI$alpha   # 7.84511
fitRI$sigma   # 11.65586
fitRI$theta   # c(3.244364e-10,  5.900000e+00, 4.300000e+00, 1.600000e+00, 2.500000e+00, 5.400000e+00)

fitRM$alpha   # 7.84511
fitRM$sigma   # 11.65586
fitRM$theta   # c(5.9, 4.3, 1.6, 2.5, 5.4)

#### Values to compare with Java classe# Pau's functions output
##### Pau's functions output
regimes <- c(1, 1, 0, 2, 2, 3, 4, 0, 0, 0, 0, 1, 2); regimes <- factor(regimes)

ouCovFI <- varOU(tr3, as.numeric(fitFI$alpha), T)
ouWFI <- weightMat(tr3, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)

ouCovFM <- varOU(tr3, as.numeric(fitFM$alpha), T)
ouWFM <- weightMat(tr3, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)

ouCovRI <- varOU(tr3, as.numeric(fitRI$alpha), F)
ouWRI <- weightMat(tr3, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)

ouCovRM <- varOU(tr3, as.numeric(fitRM$alpha), F)
ouWRM <- weightMat(tr3, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

likFI <- computeLk(dat, ouWFI %*% fitFI$theta, ouCovFI, as.numeric(fitFI$sigma))
likFM <- computeLk(dat, ouWFM %*% fitFM$theta, ouCovFM, as.numeric(fitFM$sigma))
likRI <- computeLk(dat, ouWRI %*% fitRI$theta, ouCovRI, as.numeric(fitRI$sigma))
likRM <- computeLk(dat, ouWRM %*% fitRM$theta, ouCovRM, as.numeric(fitRM$sigma))

##### LogLikelihoods
fitFI$LogLik  # -8.892192
fitFM$LogLik  # -8.892192
fitRI$LogLik  # -8.89219
fitRM$LogLik  # -8.89219

##### Comparing mvMORPH output with Pau's R functions
all.equal(as.numeric(log(likFI)), as.numeric(fitFI$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likFM)), as.numeric(fitFM$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likRI)), as.numeric(fitRI$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likRM)), as.numeric(fitRM$LogLik), tolerance = EPSILON)
