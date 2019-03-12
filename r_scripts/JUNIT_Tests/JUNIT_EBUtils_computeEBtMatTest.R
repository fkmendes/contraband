# ---------- #
# JUnit EBUtils tests expected values
# ---------- #

library(mvMORPH)

## The following three commands are aimed to modify mvMORPH package in order for it to print
## the final covariance matrix of the EB model to compare our implementation with mvMORPH's.
## We just have to compare the values of the last covariance matrix that it's shown in the display
envirMORPH <- environment(mvMORPH::mvEB)
source(file = "/Users/entimos/GitHub/contraband/r_scripts/mvEB.r")
environment(mvEB) <- envirMORPH

# The following command loads prior functions to calculate the covariance matrix,
# the weight matrix and the likelihood of the Hansen model to compare their results with mvMORPH's results
source(file = "/Users/entimos/GitHub/contraband/r_scripts/EBfunctions.R")

EPSILON <- 1e-4

# Remarks: In this script we compare both covariance matrix results and likelihood values between
# the ones provided by functions onEBfunctions.R script and the ones in mvMORPH package

# ----- START: EBUtils.computeEBtMat validation ----- #
# JUnit: EBVcvMatTest

# 1, 2, 3 refers to Test 1, Test 2 and Test 3 respectively

### Test 1: Ultrametric Tree (3 regimes)
ultTreeStr <- "(((sp1:1.0, sp2:1.0):1.0, sp3:2.0):1.0, sp4:3.0);"
ultTree <- read.tree(text = ultTreeStr)
data <- matrix(data = c(-1.3609357, -1.249852,
                        0.4767486, 2.591537,
                        0.4741497, 1.411926,
                        -0.5897473, -3.974410), nrow = 4, ncol = 2, byrow = T)
data <- data.frame(data)

  # mvEB output (mvMORPH package)
fit <- mvEB(ultTree, data, method = "inverse", param=list(up=2, low=-2))

  # Values to use on JUNIT Tests
fit$theta # c(-0.2417452, -0.3869768)
fit$beta  # 1.999999
fit$sigma
# 0.003311920 0.007755767
# 0.007755767 0.033005583

  # Values to compare with Java classes
covEB11 <- 0.66640428; covEB12 <- 0.08875625; covEB25 <- 0.2078471; # Hardcoded from mvMORPH printing
fit$LogLik  # -12.70357

  # Pau's functions output
covEB <- covmatEB(ultTree, fit$sigma, as.numeric(fit$beta))
lik1 <- computeMultiLk(data, fit$theta, covEB)

  # Comparing mvMORPH output with PAU's R functions
all.equal(as.numeric(covEB11), as.numeric(covEB[1,1]), EPSILON) # sp1 t1 vs sp1 t1
all.equal(as.numeric(covEB12), as.numeric(covEB[1,3]), EPSILON) # sp1 t1 vs sp2 t1
all.equal(as.numeric(covEB25), as.numeric(covEB[3,2]), EPSILON) # sp2 t1 vs sp1 t2

all.equal(as.numeric(log(lik1)), as.numeric(fit$LogLik), tolerance = EPSILON)


### Test 2: Non-ultrametric tree 

  # Tree data
nonultTreeStr <- "(((sp1:2.0, sp2:1.0):1.0, sp3:4.0):1.0, sp4:3.0);"
nonultTree <- read.tree(text = nonultTreeStr)
data <- matrix(data = c(-1.3609357, -1.249852,
                        0.4767486, 2.591537,
                        0.4741497, 1.411926,
                        -0.5897473, -3.974410), nrow = 4, ncol = 2)
data <- data.frame(data)

  # mvEB output (mvMORPH package)
fit <- mvEB(nonultTree, data, method = "inverse", param=list(up=2, low=-2))

  # Values to use on JUNIT Tests
fit$theta # c(0.5941884, -1.261628)
fit$beta  # 0.2204587
fit$sigma
# 0.4782103 -0.6159709
# -0.6159709  0.8544459

  # Values to compare with Java classes
covEB11 <- 3.0700716; covEB12 <- 1.2020018; covEB25 <- -1.5482690; # Hardcoded from mvMORPH printing
fit$LogLik  # -10.9459

  # Pau's functions output
covEB <- covmatEB(nonultTree, fit$sigma, as.numeric(fit$beta))
lik2 <- computeMultiLk(data, fit$theta, covEB)

  # Comparing mvMORPH output with PAU's R functions
all.equal(as.numeric(covEB11), as.numeric(covEB[1,1]), EPSILON) # sp1 t1 vs sp1 t1
all.equal(as.numeric(covEB12), as.numeric(covEB[1,3]), EPSILON) # sp1 t1 vs sp2 t1
all.equal(as.numeric(covEB25), as.numeric(covEB[3,2]), EPSILON) # sp2 t1 vs sp1 t2

all.equal(as.numeric(log(lik2)), as.numeric(fit$LogLik), tolerance = EPSILON)

# ----- END: EBUtils.computeEBtMat validation ----- #

