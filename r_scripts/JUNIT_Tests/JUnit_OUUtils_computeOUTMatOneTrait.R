# ---------- #
# JUnit OUUtils tests expected values
# ---------- #

devtools::load_all("/Users/entimos/Desktop/mvMORPH") # Modified mvMORPH package with printing of the final covariance matrix and weight matrix
source(file = "/Users/entimos/Desktop/Pau_scripts_&_functions/OUEBproject/OUfunctions.R")
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





  # Weight matrix elements (not an output of mvOU, just printed here)
wMatI12 <- 7.815427e-28; wMatI13 <- 1; wMatI32 <- 7.815427e-28; wMatI41 <- 2.184887e-41;
  # LogLikelihood
fit$LogLik  # 2.148292 

### PAU's functions output
ouCov1F <- varOU(ultTree, as.numeric(fit$alpha), T)

regimes <- c(1, 1, 2, 0, 0, 0, 1); regimes <- factor(regimes);
wMat1I <- weightMat(ultTree, as.numeric(fit$alpha), length(levels(regimes)), regimes, F)

logLikelihood1FI <- log(computeLk(ult.data, wMat1I %*% fit$theta, ouCov1F, fit$sigma))

# Comparing mvMORPH output with PAU's R functions
all.equal(ouCov1F[1,2], covMatF12, tolerance = EPSILON)
all.equal(ouCov1F[3,3], covMatF33, tolerance = EPSILON)
all.equal(ouCov1F[2,3], covMatF23, tolerance = EPSILON)
all.equal(ouCov1F[2,1], covMatF21, tolerance = EPSILON)

all.equal(wMat1I[1,2], wMatI12, tolerance = EPSILON) # Don't match accurately
all.equal(wMat1I[1,3], wMatI13, tolerance = EPSILON)
all.equal(wMat1I[3,2], wMatI32, tolerance = EPSILON)
all.equal(wMat1I[4,1], wMatI41, tolerance = EPSILON)

all.equal(as.numeric(logLikelihood1FI), as.numeric(fit$LogLik), EPSILON)


  ### FixedRoot and Merge Root case

# mvMORPH output
fit <- mvOU(ultTree, ult.data, model = "OUM", param = list(root = F, vcv = "fixedRoot"))

# Values to use on JUNIT Tests
fit$alpha   # 31.20814  
fit$sigma   # 1.248328  
fit$theta   
# group0 8.128044e-27
# group1 4.300000e+00
# group2 5.900000e+00

# Values to compare with Java classes

# Covariance matrix elements (not an output of mvOU, just printed here)
covMatF12 <- 1.563088e-29; covMatF33 <- 2.000003e-02; covMatF23 <- 1.221620e-56; covMatF21 <- 1.563088e-29;
# Weight matrix elements (not an output of mvOU, just printed here)
wMatM12 <- 1; wMatM13 <- 0; wMatM32 <- 0; wMatM41 <- 1;
# LogLikelihood
fit$LogLik  # 2.148292 

### PAU's functions output
ouCov1F <- varOU(ultTree, as.numeric(fit$alpha), T)

regimes <- c(1, 1, 2, 0, 0, 0, 1); regimes <- factor(regimes);
wMat1M <- weightMat(ultTree, as.numeric(fit$alpha), length(levels(regimes)), regimes, T)

logLikelihood1FM <- log(computeLk(ult.data, wMat1M %*% fit$theta, ouCov1F, fit$sigma))

# Comparing mvMORPH output with PAU's R functions
all.equal(ouCov1F[1,2], covMatF12, tolerance = EPSILON)
all.equal(ouCov1F[3,3], covMatF33, tolerance = EPSILON)
all.equal(ouCov1F[2,3], covMatF23, tolerance = EPSILON)
all.equal(ouCov1F[2,1], covMatF21, tolerance = EPSILON)

all.equal(wMat1M[1,2], wMatM12, tolerance = EPSILON)
all.equal(wMat1M[1,3], wMatM13, tolerance = EPSILON)
all.equal(wMat1M[3,2], wMatM32, tolerance = EPSILON)
all.equal(wMat1M[4,1], wMatM41, tolerance = EPSILON)

all.equal(as.numeric(logLikelihood1FM), as.numeric(fit$LogLik), EPSILON)


### RandomRoot and Isolated Root case

# mvMORPH output
fit <- mvOU(ultTree, ult.data, model = "OUM", param = list(root = T, vcv = "randomRoot"))

# Values to use on JUNIT Tests
fit$alpha   # 43.35287 
fit$sigma   # 1.734117 
fit$theta   
# theta_0  3.348633e-56
# group0  -1.903330e-16
# group1   4.300000e+00
# group2   5.900000e+00

# Values to compare with Java classes

# Covariance matrix elements (not an output of mvOU, just printed here)
covMatR12 <- 4.417827e-40; covMatR33 <- 2.000003e-02; covMatR23 <- 9.758589e-78; covMatR21 <- 4.417827e-40;
# Weight matrix elements (not an output of mvOU, just printed here)
wMatI12 <- 2.208911e-38; wMatI13 <- 1; wMatI32 <- 0; wMatI41 <- 1;
# LogLikelihood
fit$LogLik  # 2.148292 

### PAU's functions output
ouCov1F <- varOU(ultTree, as.numeric(fit$alpha), T)

regimes <- c(1, 1, 2, 0, 0, 0, 1); regimes <- factor(regimes);
wMat1M <- weightMat(ultTree, as.numeric(fit$alpha), length(levels(regimes)), regimes, T)

logLikelihood1FM <- log(computeLk(ult.data, wMat1M %*% fit$theta, ouCov1F, fit$sigma))

# Comparing mvMORPH output with PAU's R functions
all.equal(ouCov1F[1,2], covMatF12, tolerance = EPSILON)
all.equal(ouCov1F[3,3], covMatF33, tolerance = EPSILON)
all.equal(ouCov1F[2,3], covMatF23, tolerance = EPSILON)
all.equal(ouCov1F[2,1], covMatF21, tolerance = EPSILON)

all.equal(wMat1M[1,2], wMatM12, tolerance = EPSILON)
all.equal(wMat1M[1,3], wMatM13, tolerance = EPSILON)
all.equal(wMat1M[3,2], wMatM32, tolerance = EPSILON)
all.equal(wMat1M[4,1], wMatM41, tolerance = EPSILON)

all.equal(as.numeric(logLikelihood1FM), as.numeric(fit$LogLik), EPSILON)








  # Test 2: Non-ultrametric tree (1 regime)
nonultTreeStr <- "(((sp1:2.0, sp2:1.0):1.0, sp3:4.0):1.0, sp4:3.0);"
nonultTree <- read.tree(text = nonultTreeStr)
phylo2Mat <- vcv(nonultTree)

ouCov2F <- varOU(nonultTree, alpha, T)
ouCov2R <- varOU(nonultTree, alpha, F)

ouCov2F[1,2] # 0.08214989
ouCov2F[3,3] # 0.7136344
ouCov2F[2,3] # 0.008069795
ouCov2F[2,1] # 0.08214989

ouCov2R[1,2] # 0.08746888
ouCov2R[3,3] # 0.7142857
ouCov2R[2,3] # 0.01071113
ouCov2R[2,1] # 0.08746888

  # Test 3: Non-ultrametric tree (5 regime)
nonultTreeStrBig = "(((((sp1:1.0, sp2:1.0):1.0, sp3:1.0):2.0, (sp4:1.0, sp5:1.0):3.0):2.0, sp6:6.0):1.0, sp7:3.0);"
nonultTreeBig <- read.tree(text = nonultTreeStrBig)
phylo3Mat <- vcv(nonultTreeBig)

ouCov3F <- varOU(nonultTreeBig, alpha, T)
ouCov3R <- varOU(nonultTreeBig, alpha, F)

ouCov3F[1,2] # 0.1761011
ouCov3F[3,3] # 0.7141251
ouCov3F[2,3] # 0.08738912
ouCov3F[2,1] # 0.1761011

ouCov3R[1,2] # 0.1761407
ouCov3R[3,3] # 0.7142857
ouCov3R[2,3] # 0.08746888
ouCov3R[2,1] # 0.1761407

# ----- END: OUUtils.computeOUTMatOneTrait validation ----- #


# ----- START: OUUtils.computeWMatOneTrait validation ----- #
# JUnit: OUOneTraitcomputeWMatOneTraitTest

  # For the weight matrices:
    # 'I' suffix refers to isolating the root optimum value in the weight Matrix
    # 'M' suffix refers to merging the root parameter with the optimum parameter associated with the eldest selective regime

  # We will use the same test as before (OUOneTraitcomputeOUTMatOneTraitTest)

  # Test 1: Ultrametric Tree (3 regimes)
ultTreeStr <- "(((sp1:1.0, sp2:1.0):1.0, sp3:2.0):1.0, sp4:3.0);"
ultTree <- read.tree(text = ultTreeStr)
regimes <- c(1, 1, 2, 0, 0, 0, 1); regimes <- factor(regimes)

wMat1I <- weightMat(ultTree, alpha, length(levels(regimes)), regimes, F)
wMat1M <- weightMat(ultTree, alpha, length(levels(regimes)), regimes, T)

wMat1I[1,3] # 0.753403
wMat1I[2,2] # 0.1241405
wMat1I[3,1] # 0.1224564
wMat1I[3,3] # 0

wMat1M[1,3] # 0
wMat1M[2,2] # 0.753403
wMat1M[3,1] # 0.246597
wMat1M[3,3] # 0.753403


  # Test 2: Non-ultrametric tree (1 regime)
nonultTreeStr <- "(((sp1:2.0, sp2:1.0):1.0, sp3:4.0):1.0, sp4:3.0);"
nonultTree <- read.tree(text = nonultTreeStr)
regimes <- c(0, 0, 0, 0, 0, 0, 0); regimes <- factor(regimes)

wMat2I <- weightMat(nonultTree, alpha, length(levels(regimes)), regimes, F)
wMat2M <- weightMat(nonultTree, alpha, length(levels(regimes)), regimes, T)

wMat2I[1,1] # 0.06081006
wMat2I[1,2] # 0.9391899
wMat2I[3,1] # 0.03019738
wMat2I[3,2] # 0.9698026

wMat2M[1,1] # 1
wMat2M[2,1] # 1
wMat2M[3,1] # 1
wMat2M[4,1] # 1

# Test 3: Non-ultrametric tree (5 regime)
nonultTreeStrBig <- "(((((sp1:1.0, sp2:1.0):1.0, sp3:1.0):2.0, (sp4:1.0, sp5:1.0):3.0):2.0, sp6:6.0):1.0, sp7:3.0);"
nonultTreeBig <- read.tree(text = nonultTreeStrBig)
regimes <- c(1, 1, 0, 2, 2, 3, 4, 0, 0, 0, 0, 1, 2); regimes <- factor(regimes)

wMat3I <- weightMat(nonultTreeBig, alpha, length(levels(regimes)), regimes, F)
wMat3M <- weightMat(nonultTreeBig, alpha, length(levels(regimes)), regimes, T)

wMat3I[1,1] # 0.007446583
wMat3I[1,2] # 0.239150381
wMat3I[3,1] # 0.014995577
wMat3I[3,2] # 0.985004423

wMat3M[1,1] # 0.24659696
wMat3M[2,1] # 0.24659696
wMat3M[3,1] # 1.00000000
wMat3M[4,1] # 0.06081006

# ----- END: OUUtils.computeWMatOneTrait validation ----- #

