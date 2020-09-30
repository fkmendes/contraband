# author: Pau Bravo
# This R script gives us the expected values for JUnit tests
# 'OUUtilsComputeWMatOneTraitTest'
# (1) 'testComputeOUWMatSmallTree'
# (2) 'testComputeOUWMatSmallTreeNonUltra'
# (3) 'testComputeOUWMatLargerTreeNonUltra'

# devtools::load_all("/path/to/modified/mvMORPH") # Modified mvMORPH package with printing of the final covariance matrix and weight matrix

## The following three commands are aimed to modify mvMORPH package in order for it to print
## the final covariance matrix (multiplied by sigma^2) and weight matrix of the OU model to compare our implementation with mvMORPH's.
envirMORPH <- environment(mvMORPH::mvOU)

# The following command loads prior functions to calculate the covariance matrix,
# the weight matrix and the likelihood of the Hansen model to compare their results with mvMORPH's results
source(file="../OUfunctions.R")

EPSILON <- 1e-04 # Tolerance value for comparing mvMORPH with PAU's functions

## 'OUUtilsComputeWMatOneTraitTest'

## R: rootIsRandVar=true (vcv="randomRoot")
## F: rootIsRandVar=false (vcv="fixedRoot")
## I: useRootMetaData=true (root=T)
## M: useRootMetaData=false (root=F)

### (1) testComputeOUWMatSmallTree

tr <- paintSubTree(read.tree(text="(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);"), node=6, state="group1", anc.state="group0")
tr <- paintSubTree(tr, node=3, state="group2", anc.state="group0", stem=T) # note that in the Java tests, we do not need to care about regimes (we will be testing OUUtils given all the MLEs)
ult.data <- c(4.1, 4.5, 5.9, 0.0)
plotSimmap(tr)

fitFI <- mvOU(tr, ult.data, model="OUM", param=list(root=T, vcv="fixedRoot"))
fitFM <- mvOU(tr, ult.data, model="OUM", param=list(root=F, vcv="fixedRoot"))
fitRI <- mvOU(tr, ult.data, model="OUM", param=list(root=T, vcv="randomRoot"))
fitRM <- mvOU(tr, ult.data, model="OUM", param=list(root=F, vcv="randomRoot"))

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

#### Values to compare with Java classes
##### Pau's functions output
regimes <- c(1, 1, 2, 0, 0, 0, 1); regimes <- factor(regimes);
ouWFI <- weightMat(tr, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)
ouWFM <- weightMat(tr, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)
ouWRI <- weightMat(tr, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)
ouWRM <- weightMat(tr, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

##### Hardcoded expected weight matrix elements (not an output of mvOU, just printed here)
wMatFI12 <- 7.815427e-28; wMatFI33 <- 0; wMatFI23 <- 1; wMatFI21 <- 2.184887e-41;
wMatFM12 <- 1;            wMatFM33 <- 1; wMatFM23 <- 0; wMatFM21 <- 7.815427e-28;
wMatRI12 <- 2.208911e-38; wMatRI33 <- 0; wMatRI23 <- 1; wMatRI21 <- 3.282974e-57;
wMatRM12 <- 1;            wMatRM33 <- 1; wMatRM23 <- 0; wMatRM21 <- 2.208911e-38;

##### Comparing mvMORPH output with PAU's R functions
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

### (2) testComputeOUWMatSmallTreeNonUltra

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

#### Values to compare with Java classes
##### Pau's functions output
regimes <- c(0, 0, 0, 0, 0, 0, 0); regimes <- factor(regimes);
ouWFI <- weightMat(tr2, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)
ouWFM <- weightMat(tr2, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)
ouWRI <- weightMat(tr2, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)
ouWRM <- weightMat(tr2, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

##### Hardcoded expected weight matrix elements (not an output of mvOU, just printed here)
wMatFI12 <- 0.9495264;  wMatFI22 <- 0.8935126;  wMatFI31 <- 0.02392380; wMatFI41 <- 0.10648738;
wMatFM11 <- 1;          wMatFM21 <- 1;          wMatFM31 <- 1;          wMatFM41 <- 1;
wMatRI12 <- 0.9680612;  wMatRI22 <- 0.9244492;  wMatRI31 <- 0.01350201; wMatRI41 <- 0.07555081;
wMatRM11 <- 1;          wMatRM21 <- 1;          wMatRM31 <- 1;          wMatRM41 <- 1;

##### Comparing mvMORPH output with Pau's R functions
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

### (3) testComputeOUWMatLargerTreeNonUltra

tr3 <- paintSubTree(read.tree(text="(((((sp1:1.0, sp2:1.0):1.0, sp3:1.0):2.0, (sp4:1.0, sp5:1.0):3.0):2.0, sp6:6.0):1.0, sp7:3.0);"),
                              node=12, state="group1", anc.state="group0", stem=T)
tr3 <- paintSubTree(tr3, node=12, state="group1", anc.state="group0", stem=T)
tr3 <- paintSubTree(tr3, node=13, state="group2", anc.state="group0", stem=T)
tr3 <- paintSubTree(tr3, node=6, state="group3", anc.state="group0", stem=T)
tr3 <- paintSubTree(tr3, node=7, state="group4", anc.state="group0", stem=T)
dat2 <- c(4.1, 4.5, 5.9, 0.0, 3.2, 2.5, 5.4)
plotSimmap(tr3)

fitFI <- mvOU(tr3, dat2, model = "OUM", param = list(root = T, vcv = "fixedRoot"))
fitFM <- mvOU(tr3, dat2, model = "OUM", param = list(root = F, vcv = "fixedRoot"))
fitRI <- mvOU(tr3, dat2, model = "OUM", param = list(root = T, vcv = "randomRoot"))
fitRM <- mvOU(tr3, dat2, model = "OUM", param = list(root = F, vcv = "randomRoot"))

#### Values to use on JUNIT Tests
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

#### Values to compare with Java classes
##### Pau's functions output
regimes <- c(1, 1, 0, 2, 2, 3, 4, 0, 0, 0, 0, 1, 2); regimes <- factor(regimes)
ouWFI <- weightMat(tr3, as.numeric(fitFI$alpha), length(levels(regimes)), regimes, F)
ouWFM <- weightMat(tr3, as.numeric(fitFM$alpha), length(levels(regimes)), regimes, T)
ouWRI <- weightMat(tr3, as.numeric(fitRI$alpha), length(levels(regimes)), regimes, F)
ouWRM <- weightMat(tr3, as.numeric(fitRM$alpha), length(levels(regimes)), regimes, T)

##### Hardcoded expected weight matrix elements (not an output of mvOU, just printed here)
wMatFI12 <- 6.247136e-07; wMatFI33 <- 0; wMatFI23 <- 0.9999994; wMatFI21 <- 1.927008e-22;
wMatFM12 <- 0.9999994;    wMatFM33 <- 0; wMatFM23 <- 0;         wMatFM21 <- 6.247136e-07;
wMatRI12 <- 1.533995e-07; wMatRI33 <- 0; wMatRI23 <- 0.9999998; wMatRI21 <- 1.413785e-24;
wMatRM12 <- 0.9999998;    wMatRM33 <- 0; wMatRM23 <- 0;         wMatRM21 <- 1.533995e-07;

# Comparing mvMORPH output with Pau's R functions
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

