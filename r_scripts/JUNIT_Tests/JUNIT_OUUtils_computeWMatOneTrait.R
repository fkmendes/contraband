# ---------- #
# JUnit OUUtils tests expected values
# ---------- #

devtools::load_all("/Users/entimos/Desktop/mvMORPH") # Modified mvMORPH package with printing of the final covariance matrix and weight matrix
source(file = "/Users/entimos/Desktop/Pau_scripts_&_functions/OUEBproject/OUfunctions.R")
EPSILON <- 1e-04 # Tolerance value for comparing mvMORPH with PAU's functions

# 1, 2, 3 refers to Test 1, Test 2 and Test 3 respectively

# ----- START: OUUtils.computeWMatOneTrait validation ----- #
# JUnit: OUOneTraitcomputeWMatOneTraitTest 

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

# Values to compare with Java classes

# weight matrix elements (not an output of mvOU, just printed here)
wMatFI12 <- 7.815427e-28; wMatFI33 <- 0; wMatFI23 <- 1; wMatFI21 <- 2.184887e-41;
wMatFM12 <- 1;            wMatFM33 <- 1; wMatFM23 <- 0; wMatFM21 <- 7.815427e-28;
wMatRI12 <- 2.208911e-38; wMatRI33 <- 0; wMatRI23 <- 1; wMatRI21 <- 3.282974e-57;
wMatRM12 <- 1;            wMatRM33 <- 1; wMatRM23 <- 0; wMatRM21 <- 2.208911e-38;

# PAU's functions output
regimes <- c(1, 1, 2, 0, 0, 0, 1); regimes <- factor(regimes);
ouWFI <- weightMat(ultTree, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)
ouWFM <- weightMat(ultTree, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)
ouWRI <- weightMat(ultTree, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)
ouWRM <- weightMat(ultTree, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

# Comparing mvMORPH output with PAU's R functions
all.equal(ouWFI[1,2], wMatFI12, tolerance = EPSILON)
all.equal(ouWFI[3,3], wMatFI33, tolerance = EPSILON) 
all.equal(ouWFI[2,3], wMatFI23, tolerance = EPSILON)
all.equal(ouWFI[2,1], wMatFI21, tolerance = EPSILON)

all.equal(ouWFM[1,2], wMatFM12, tolerance = EPSILON)
all.equal(ouWFM[3,3], wMatFM33, tolerance = EPSILON)
all.equal(ouWFM[2,3], wMatFM23, tolerance = EPSILON)
all.equal(ouWFM[2,1], wMatFM21, tolerance = EPSILON)

all.equal(ouWRI[1,2], wMatRI12, tolerance = EPSILON)
all.equal(ouWRI[3,3], wMatRI33, tolerance = EPSILON)
all.equal(ouWRI[2,3], wMatRI23, tolerance = EPSILON)
all.equal(ouWRI[2,1], wMatRI21, tolerance = EPSILON)

all.equal(ouWRM[1,2], wMatRM12, tolerance = EPSILON)
all.equal(ouWRM[3,3], wMatRM33, tolerance = EPSILON)
all.equal(ouWRM[2,3], wMatRM23, tolerance = EPSILON)
all.equal(ouWRM[2,1], wMatRM21, tolerance = EPSILON)


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
wMatFI12 <- 0.9495264;  wMatFI22 <- 0.8935126;  wMatFI31 <- 0.02392380; wMatFI41 <- 0.10648738;
wMatFM11 <- 1;          wMatFM21 <- 1;          wMatFM31 <- 1;          wMatFM41 <- 1;
wMatRI12 <- 0.9680612;  wMatRI22 <- 0.9244492;  wMatRI31 <- 0.01350201; wMatRI41 <- 0.07555081;
wMatRM11 <- 1;          wMatRM21 <- 1;          wMatRM31 <- 1;          wMatRM41 <- 1;

# PAU's functions output
regimes <- c(0, 0, 0, 0, 0, 0, 0); regimes <- factor(regimes);
ouWFI <- weightMat(nonultTree, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)
ouWFM <- weightMat(nonultTree, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)
ouWRI <- weightMat(nonultTree, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)
ouWRM <- weightMat(nonultTree, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

# Comparing mvMORPH output with PAU's R functions
all.equal(ouWFI[1,2], wMatFI12, tolerance = EPSILON)
all.equal(ouWFI[2,2], wMatFI22, tolerance = EPSILON)  
all.equal(ouWFI[3,1], wMatFI31, tolerance = EPSILON)
all.equal(ouWFI[4,1], wMatFI41, tolerance = EPSILON)

all.equal(ouWFM[1,1], wMatFM11, tolerance = EPSILON)
all.equal(ouWFM[2,1], wMatFM21, tolerance = EPSILON)  
all.equal(ouWFM[3,1], wMatFM31, tolerance = EPSILON)
all.equal(ouWFM[4,1], wMatFM41, tolerance = EPSILON)

all.equal(ouWRI[1,2], wMatRI12, tolerance = EPSILON)
all.equal(ouWRI[2,2], wMatRI22, tolerance = EPSILON)  
all.equal(ouWRI[3,1], wMatRI31, tolerance = EPSILON)
all.equal(ouWRI[4,1], wMatRI41, tolerance = EPSILON)

all.equal(ouWRM[1,1], wMatRM11, tolerance = EPSILON)
all.equal(ouWRM[2,1], wMatRM21, tolerance = EPSILON)  
all.equal(ouWRM[3,1], wMatRM31, tolerance = EPSILON)
all.equal(ouWRM[4,1], wMatRM41, tolerance = EPSILON)


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
wMatFI12 <- 6.247136e-07; wMatFI33 <- 0; wMatFI23 <- 0.9999994; wMatFI21 <- 1.927008e-22;
wMatFM12 <- 0.9999994;    wMatFM33 <- 0; wMatFM23 <- 0;         wMatFM21 <- 6.247136e-07;
wMatRI12 <- 1.533995e-07; wMatRI33 <- 0; wMatRI23 <- 0.9999998; wMatRI21 <- 1.413785e-24;
wMatRM12 <- 0.9999998;    wMatRM33 <- 0; wMatRM23 <- 0;         wMatRM21 <- 1.533995e-07;

# PAU's functions output
regimes <- c(1, 1, 0, 2, 2, 3, 4, 0, 0, 0, 0, 1, 2); regimes <- factor(regimes)
ouWFI <- weightMat(nonultTree, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)
ouWFM <- weightMat(nonultTree, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)
ouWRI <- weightMat(nonultTree, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)
ouWRM <- weightMat(nonultTree, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

# Comparing mvMORPH output with PAU's R functions
all.equal(ouWFI[1,2], wMatFI12, tolerance = EPSILON)
all.equal(ouWFI[3,3], wMatFI33, tolerance = EPSILON)  
all.equal(ouWFI[2,3], wMatFI23, tolerance = EPSILON) # "Mean relative difference: 0.000790413"
all.equal(ouWFI[2,1], wMatFI21, tolerance = EPSILON)

all.equal(ouWFM[1,2], wMatFM12, tolerance = EPSILON)
all.equal(ouWFM[3,3], wMatFM33, tolerance = EPSILON)  
all.equal(ouWFM[2,3], wMatFM23, tolerance = EPSILON)
all.equal(ouWFM[2,1], wMatFM21, tolerance = EPSILON)

all.equal(ouWRI[1,2], wMatRI12, tolerance = EPSILON)
all.equal(ouWRI[3,3], wMatRI33, tolerance = EPSILON)  
all.equal(ouWRI[2,3], wMatRI23, tolerance = EPSILON) # "Mean relative difference: 0.0003916158"
all.equal(ouWRI[2,1], wMatRI21, tolerance = EPSILON)

all.equal(ouWRM[1,2], wMatRM12, tolerance = EPSILON)
all.equal(ouWRM[3,3], wMatRM33, tolerance = EPSILON)  
all.equal(ouWRM[2,3], wMatRM23, tolerance = EPSILON)
all.equal(ouWRM[2,1], wMatRM21, tolerance = EPSILON)

# ----- END: OUUtils.computeWMatOneTrait validation ----- #

