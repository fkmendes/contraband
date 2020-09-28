# ---------- #
# JUnit OUUtils tests expected values
# ---------- #

library(mvMORPH)

## The following three commands are aimed to modify mvMORPH package in order for it to print
## the final covariance matrix (multiplied by sigma^2) and weight matrix of the OU model to compare our implementation with mvMORPH's.
envirMORPH <- environment(mvMORPH::mvOU)
source(file = "/Users/entimos/GitHub/contraband/r_scripts/mvOU.r")
environment(mvOU) <- envirMORPH

# The following command loads prior functions to calculate the covariance matrix,
# the weight matrix and the likelihood of the Hansen model to compare their results with mvMORPH's results
source(file = "/Users/entimos/GitHub/contraband/r_scripts/OUfunctions.R")

EPSILON <- 1e-04 # Tolerance value for comparing mvMORPH with PAU's functions

# 1, 2, 3 refers to Test 1, Test 2 and Test 3 respectively

# ----- START: MVNUtils.computeMVNLk validation ----- #
# JUnit: OUOneTraitcomputeMVNLkTest 

# For the covariance matrices we can assume two different hipothesis:
# 'F' suffix refers to assuming a fixed parameter Root
# 'R' suffix refers to assuming a random variable Root (stationary distribution)
# For the weight matrices we can assume two different hipothesis:
# 'I' suffix refers to isolating the root optimum value in the weight Matrix
# 'M' suffix refers to merging the root parameter with the optimum parameter associated with the eldest selective regime

# Every test has four different outputs according to the combination of the previous situations: FI, FM, RI, RM

### Test 1: Ultrametric Tree (3 regimes)

# Tree data
ultTreeStr <- "(((sp1:1.0, sp2:1.0):1.0, sp3:2.0):1.0, sp4:3.0);"
ultTree <- read.tree(text = ultTreeStr)
ultTree <- paintSubTree(ultTree, node = 6, state = "group1", anc.state = "group0")
ultTree <- paintSubTree(ultTree, node = 3, state = "group2", anc.state = "group0", stem = T)
ult.data <- c(4.1, 4.5, 5.9, 0.0)
plot(ultTree)

# mvMORPH output
fitFI <- mvOU(ultTree, ult.data, model = "OUM", param = list(root = T, vcv = "fixedRoot"))
fitFM <- mvOU(ultTree, ult.data, model = "OUM", param = list(root = F, vcv = "fixedRoot"))
fitRI <- mvOU(ultTree, ult.data, model = "OUM", param = list(root = T, vcv = "randomRoot"))
fitRM <- mvOU(ultTree, ult.data, model = "OUM", param = list(root = F, vcv = "randomRoot"))

# Values to use on JUNIT Tests

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

# LogLikelihood values to compare with Java classes
fitFI$LogLik  # 2.148292
fitFM$LogLik  # 2.148292
fitRI$LogLik  # 2.148292
fitRM$LogLik  # 2.148292

# PAU's functions output
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

# Comparing mvMORPH output with PAU's R functions
all.equal(as.numeric(log(likFI)), as.numeric(fitFI$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likFM)), as.numeric(fitFM$LogLik), tolerance = EPSILON) 
all.equal(as.numeric(log(likRI)), as.numeric(fitRI$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likRM)), as.numeric(fitRM$LogLik), tolerance = EPSILON)


### Test 2: Non-ultrametric tree (1 regime)

# Tree data
nonultTreeStr <- "(((sp1:2.0, sp2:1.0):1.0, sp3:4.0):1.0, sp4:3.0);"
nonultTree <- read.tree(text = nonultTreeStr)
nonult.data <- c(4.1, 4.5, 5.9, 0.0)
plot(nonultTree)

# mvMORPH output
fitFI <- mvOU(nonultTree, nonult.data, model = "OU1", param = list(root = T, vcv = "fixedRoot"))
fitFM <- mvOU(nonultTree, nonult.data, model = "OU1", param = list(root = F, vcv = "fixedRoot"))
fitRI <- mvOU(nonultTree, nonult.data, model = "OU1", param = list(root = T, vcv = "randomRoot"))
fitRM <- mvOU(nonultTree, nonult.data, model = "OU1", param = list(root = F, vcv = "randomRoot"))

# Values to use on JUNIT Tests

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

# LogLikelihood values to compare with Java classes
fitFI$LogLik  # -7.630117
fitFM$LogLik  # -8.457486
fitRI$LogLik  # -7.63854
fitRM$LogLik  # -8.817273

# PAU's functions output
regimes <- c(0, 0, 0, 0, 0, 0, 0); regimes <- factor(regimes);

ouCovFI <- varOU(nonultTree, as.numeric(fitFI$alpha), T)
ouWFI <- weightMat(nonultTree, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)

ouCovFM <- varOU(nonultTree, as.numeric(fitFM$alpha), T)
ouWFM <- weightMat(nonultTree, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)

ouCovRI <- varOU(nonultTree, as.numeric(fitRI$alpha), F)
ouWRI <- weightMat(nonultTree, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)

ouCovRM <- varOU(nonultTree, as.numeric(fitRM$alpha), F)
ouWRM <- weightMat(nonultTree, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

likFI <- computeLk(nonult.data, ouWFI %*% fitFI$theta, ouCovFI, as.numeric(fitFI$sigma))
likFM <- computeLk(nonult.data, ouWFM %*% fitFM$theta, ouCovFM, as.numeric(fitFM$sigma))
likRI <- computeLk(nonult.data, ouWRI %*% fitRI$theta, ouCovRI, as.numeric(fitRI$sigma))
likRM <- computeLk(nonult.data, ouWRM %*% fitRM$theta, ouCovRM, as.numeric(fitRM$sigma))

# Comparing mvMORPH output with PAU's R functions
all.equal(as.numeric(log(likFI)), as.numeric(fitFI$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likFM)), as.numeric(fitFM$LogLik), tolerance = EPSILON) 
all.equal(as.numeric(log(likRI)), as.numeric(fitRI$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likRM)), as.numeric(fitRM$LogLik), tolerance = EPSILON)


### Test 3: Non-ultrametric tree (5 regimes)

# Tree data
nonultTreeStrBig = "(((((sp1:1.0, sp2:1.0):1.0, sp3:1.0):2.0, (sp4:1.0, sp5:1.0):3.0):2.0, sp6:6.0):1.0, sp7:3.0);"
nonultTreeBig <- read.tree(text = nonultTreeStrBig)
nonultTreeBig <- paintSubTree(nonultTreeBig, node = 12, state = "group1", anc.state = "group0", stem = T)
nonultTreeBig <- paintSubTree(nonultTreeBig, node = 13, state = "group2", anc.state = "group0", stem = T)
nonultTreeBig <- paintSubTree(nonultTreeBig, node = 6, state = "group3", anc.state = "group0", stem = T)
nonultTreeBig <- paintSubTree(nonultTreeBig, node = 7, state = "group4", anc.state = "group0", stem = T)
nonultTreeBig.data <- c(4.1, 4.5, 5.9, 0.0, 3.2, 2.5, 5.4)
plot(nonultTreeBig)

# mvMORPH output
fitFI <- mvOU(nonultTreeBig, nonultTreeBig.data, model = "OUM", param = list(root = T, vcv = "fixedRoot"))
fitFM <- mvOU(nonultTreeBig, nonultTreeBig.data, model = "OUM", param = list(root = F, vcv = "fixedRoot"))
fitRI <- mvOU(nonultTreeBig, nonultTreeBig.data, model = "OUM", param = list(root = T, vcv = "randomRoot"))
fitRM <- mvOU(nonultTreeBig, nonultTreeBig.data, model = "OUM", param = list(root = F, vcv = "randomRoot"))

# Values to use on JUNIT Tests

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


# LogLikelihood values to compare with Java classes
fitFI$LogLik  # -8.892192
fitFM$LogLik  # -8.892192
fitRI$LogLik  # -8.89219
fitRM$LogLik  # -8.89219

# PAU's functions output
regimes <- c(1, 1, 0, 2, 2, 3, 4, 0, 0, 0, 0, 1, 2); regimes <- factor(regimes)

ouCovFI <- varOU(nonultTreeBig, as.numeric(fitFI$alpha), T)
ouWFI <- weightMat(nonultTreeBig, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)

ouCovFM <- varOU(nonultTreeBig, as.numeric(fitFM$alpha), T)
ouWFM <- weightMat(nonultTreeBig, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)

ouCovRI <- varOU(nonultTreeBig, as.numeric(fitRI$alpha), F)
ouWRI <- weightMat(nonultTreeBig, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)

ouCovRM <- varOU(nonultTreeBig, as.numeric(fitRM$alpha), F)
ouWRM <- weightMat(nonultTreeBig, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

likFI <- computeLk(nonultTreeBig.data, ouWFI %*% fitFI$theta, ouCovFI, as.numeric(fitFI$sigma))
likFM <- computeLk(nonultTreeBig.data, ouWFM %*% fitFM$theta, ouCovFM, as.numeric(fitFM$sigma))
likRI <- computeLk(nonultTreeBig.data, ouWRI %*% fitRI$theta, ouCovRI, as.numeric(fitRI$sigma))
likRM <- computeLk(nonultTreeBig.data, ouWRM %*% fitRM$theta, ouCovRM, as.numeric(fitRM$sigma))

# Comparing mvMORPH output with PAU's R functions
all.equal(as.numeric(log(likFI)), as.numeric(fitFI$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likFM)), as.numeric(fitFM$LogLik), tolerance = EPSILON) 
all.equal(as.numeric(log(likRI)), as.numeric(fitRI$LogLik), tolerance = EPSILON)
all.equal(as.numeric(log(likRM)), as.numeric(fitRM$LogLik), tolerance = EPSILON)

# ----- END: MVNUtils.computeMVNLk validation ----- #
