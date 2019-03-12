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

# ----- START: OUUtils.computeOUTMatOneTrait validation ----- #
# JUnit: OUOneTraitcomputeOUTMatOneTraitTest 

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
fitFM$sigma   # 1.248328 
fitFM$theta   # c(group0, group1, group2) = c(8.128044e-27, 4.300000e+00, 5.900000e+00)

fitRI$alpha   # 43.35287
fitRI$sigma   # 1.734117
fitRI$theta   # c(theta0, group0, group1, group2) = c(3.348633e-56, -1.903330e-16, 4.300000e+00, 5.900000e+00)

fitRM$alpha   # 43.35287
fitRM$sigma   # 1.734117 
fitRM$theta   # c(group0, group1, group2) = c(2.297268e-37, 4.300000e+00, 5.900000e+00)

# Values to compare with Java classes

  # Covariance matrix elements (not an output of mvOU, just printed here)
covMatFI12 <- 1.563088e-29; covMatFI33 <- 2.000003e-02; covMatFI23 <- 1.221620e-56; covMatFI21 <- 1.563088e-29;
covMatFM12 <- 1.563088e-29; covMatFM33 <- 2.000003e-02; covMatFM23 <- 1.221620e-56; covMatFM21 <- 1.563088e-29;
covMatRI12 <- 4.417827e-40; covMatRI33 <- 2.000002e-02; covMatRI23 <- 9.758589e-78; covMatRI21 <- 4.417827e-40;
covMatRM12 <- 4.417827e-40; covMatRM33 <- 2.000002e-02; covMatRM23 <- 9.758589e-78; covMatRM21 <- 4.417827e-40;

  # PAU's functions output
ouCovFI <- as.numeric(fitFI$sigma) * varOU(ultTree, as.numeric(fitFI$alpha), T)
ouCovFM <- as.numeric(fitFM$sigma) * varOU(ultTree, as.numeric(fitFM$alpha), T)
ouCovRI <- as.numeric(fitRI$sigma) * varOU(ultTree, as.numeric(fitRI$alpha), F)
ouCovRM <- as.numeric(fitRM$sigma) * varOU(ultTree, as.numeric(fitRM$alpha), F)

  # Comparing mvMORPH output with PAU's R functions
all.equal(ouCovFI[1,2], covMatFI12, tolerance = EPSILON)
all.equal(ouCovFI[3,3], covMatFI33, tolerance = EPSILON) 
all.equal(ouCovFI[2,3], covMatFI23, tolerance = EPSILON)
all.equal(ouCovFI[2,1], covMatFI21, tolerance = EPSILON)

all.equal(ouCovFM[1,2], covMatFM12, tolerance = EPSILON)
all.equal(ouCovFM[3,3], covMatFM33, tolerance = EPSILON)
all.equal(ouCovFM[2,3], covMatFM23, tolerance = EPSILON)
all.equal(ouCovFM[2,1], covMatFM21, tolerance = EPSILON)

all.equal(ouCovRI[1,2], covMatRI12, tolerance = EPSILON)
all.equal(ouCovRI[3,3], covMatRI33, tolerance = EPSILON)
all.equal(ouCovRI[2,3], covMatRI23, tolerance = EPSILON)
all.equal(ouCovRI[2,1], covMatRI21, tolerance = EPSILON)

all.equal(ouCovRM[1,2], covMatRM12, tolerance = EPSILON)
all.equal(ouCovRM[3,3], covMatRM33, tolerance = EPSILON)
all.equal(ouCovRM[2,3], covMatRM23, tolerance = EPSILON)
all.equal(ouCovRM[2,1], covMatRM21, tolerance = EPSILON)


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

# Values to compare with Java classes

# Covariance matrix elements (not an output of mvOU, just printed here)
covMatFI12 <- 0.2711105;  covMatFI33 <- 2.67973939;   covMatFI23 <- 0.02357370;   covMatFI21 <- 0.27111053;
covMatFM12 <- 2.475727;   covMatFM33 <- 6.189318;     covMatFM23 <- 1.237863;     covMatFM21 <- 2.475727;
covMatRI12 <- 0.20187480; covMatRI33 <- 2.672040045;  covMatRI23 <- 0.015251805;  covMatRI21 <- 0.201874800;
covMatRM12 <- 0.57628838; covMatRM33 <- 4.82816083;   covMatRM23 <- 0.06878567;   covMatRM21 <- 0.57628838;

# PAU's functions output
ouCovFI <- as.numeric(fitFI$sigma) * varOU(nonultTree, as.numeric(fitFI$alpha), T)
ouCovFM <- as.numeric(fitFM$sigma) * varOU(nonultTree, as.numeric(fitFM$alpha), T)
ouCovRI <- as.numeric(fitRI$sigma) * varOU(nonultTree, as.numeric(fitRI$alpha), F)
ouCovRM <- as.numeric(fitRM$sigma) * varOU(nonultTree, as.numeric(fitRM$alpha), F)

# Comparing mvMORPH output with PAU's R functions
all.equal(ouCovFI[1,2], covMatFI12, tolerance = EPSILON)
all.equal(ouCovFI[3,3], covMatFI33, tolerance = EPSILON)  
all.equal(ouCovFI[2,3], covMatFI23, tolerance = EPSILON)
all.equal(ouCovFI[2,1], covMatFI21, tolerance = EPSILON)

all.equal(ouCovFM[1,2], covMatFM12, tolerance = EPSILON)
all.equal(ouCovFM[3,3], covMatFM33, tolerance = EPSILON)  
all.equal(ouCovFM[2,3], covMatFM23, tolerance = EPSILON)
all.equal(ouCovFM[2,1], covMatFM21, tolerance = EPSILON)

all.equal(ouCovRI[1,2], covMatRI12, tolerance = EPSILON)
all.equal(ouCovRI[3,3], covMatRI33, tolerance = EPSILON)  
all.equal(ouCovRI[2,3], covMatRI23, tolerance = EPSILON)
all.equal(ouCovRI[2,1], covMatRI21, tolerance = EPSILON)

all.equal(ouCovRM[1,2], covMatRM12, tolerance = EPSILON)
all.equal(ouCovRM[3,3], covMatRM33, tolerance = EPSILON)  
all.equal(ouCovRM[2,3], covMatRM23, tolerance = EPSILON)
all.equal(ouCovRM[2,1], covMatRM21, tolerance = EPSILON)


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

# Values to compare with Java classes

# Covariance matrix elements (not an output of mvOU, just printed here)
covMatFI12 <- 4.640928e-07; covMatFI33 <- 7.428888e-01; covMatFI23 <- 3.668135e-10; covMatFI21 <- 4.640928e-07;
covMatFM12 <- 4.640928e-07; covMatFM33 <- 7.428888e-01; covMatFM23 <- 3.668135e-10; covMatFM21 <- 4.640928e-07;
covMatRI12 <- 1.139565e-07; covMatRI33 <- 7.428742e-01; covMatRI23 <- 4.463248e-11; covMatRI21 <- 1.139565e-07;
covMatRM12 <- 1.139565e-07; covMatRM33 <- 7.428742e-01; covMatRM23 <- 4.463248e-11; covMatRM21 <- 1.139565e-07;

# PAU's functions output
ouCovFI <- as.numeric(fitFI$sigma) * varOU(nonultTreeBig, as.numeric(fitFI$alpha), T)
ouCovFM <- as.numeric(fitFM$sigma) * varOU(nonultTreeBig, as.numeric(fitFM$alpha), T)
ouCovRI <- as.numeric(fitRI$sigma) * varOU(nonultTreeBig, as.numeric(fitRI$alpha), T)
ouCovRM <- as.numeric(fitRM$sigma) * varOU(nonultTreeBig, as.numeric(fitRM$alpha), T)

# Comparing mvMORPH output with PAU's R functions
all.equal(ouCovFI[1,2], covMatFI12, tolerance = EPSILON)
all.equal(ouCovFI[3,3], covMatFI33, tolerance = EPSILON)  
all.equal(ouCovFI[2,3], covMatFI23, tolerance = EPSILON)
all.equal(ouCovFI[2,1], covMatFI21, tolerance = EPSILON)

all.equal(ouCovFM[1,2], covMatFM12, tolerance = EPSILON)
all.equal(ouCovFM[3,3], covMatFM33, tolerance = EPSILON)  
all.equal(ouCovFM[2,3], covMatFM23, tolerance = EPSILON)
all.equal(ouCovFM[2,1], covMatFM21, tolerance = EPSILON)

all.equal(ouCovRI[1,2], covMatRI12, tolerance = EPSILON)
all.equal(ouCovRI[3,3], covMatRI33, tolerance = EPSILON)  
all.equal(ouCovRI[2,3], covMatRI23, tolerance = EPSILON)
all.equal(ouCovRI[2,1], covMatRI21, tolerance = EPSILON)

all.equal(ouCovRM[1,2], covMatRM12, tolerance = EPSILON)
all.equal(ouCovRM[3,3], covMatRM33, tolerance = EPSILON) 
all.equal(ouCovRM[2,3], covMatRM23, tolerance = EPSILON)
all.equal(ouCovRM[2,1], covMatRM21, tolerance = EPSILON)

# ----- START: OUUtils.computeOUTMatOneTrait validation ----- #



 