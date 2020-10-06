# This R script gives us the expected values for JUnit tests
# 'NodeMathTest'
# (1) 'testShrinkageOperations'
# (2) 'testNonShrinkageOperationsOneRateOnly'
# (3) 'testNonShrinkageOperationsMultipleRates'

library(mcmc3r)
library(corpcor)

## 'NodeMathTest'

### (1) testShrinkageOperations
dat <- matrix(c(1, 3, 2, 2, 5, 4), nrow = 3, ncol = 2)

options(digits = 16)
unBiasedEstimation <- cor(dat)[1,2] # 0.9819805060619656

# shrinkage parameter
options(digits = 15)
delta <- estimate.lambda(x = dat, verbose = FALSE)

#### shrinkage estimation
shrinkage.correlation.matrix <- matrix(cor.shrink(dat, verbose = FALSE), 2, 2)
# shrinkage.correlation.matrix[1,2] -> 0.727392967453307

#### determinant of shrinkage correlation matrix and inverse shrinkage correlation matrix
det.rho <- det(shrinkage.correlation.matrix) # 0.470899470899473
det.inv.rho <- det(solve(shrinkage.correlation.matrix))

#### trait rate matrix
sigmasq <- 0.1543038
trait.rate.matrix <- sigmasq * shrinkage.correlation.matrix
# trait.rate.matrix[1,2] -> 0.1122394989713215 

#### determinant of trait rate matrix and inverse trait rate matrix
det.trait.rate.matrix <- log(det(trait.rate.matrix)) # -4.4907744304614
det.inv.traitRate.matrix <- log(det(solve(trait.rate.matrix))) # 4.4907744304614

#### transform the data
lMat = chol(trait.rate.matrix) 
# inverse of lMat -> A = L.inverse
aMat = solve(lMat)
# transformed data -> Z = M * A
Z = dat %*% aMat # 2.54572618261879, 4.72108594804791, 7.63717854785638, 10.4534826917935, 5.09145236523758, 9.44217189609582

### (2) testNonShrinkageOperationsOneRateOnly
dat <- matrix(c(-2.62762948691895, -1.56292164859448,
                -1.50846427625826, -1.59482814741543,
                -0.226074849617958, -2.11000367246907), nrow = 3, ncol = 2, byrow = TRUE)

#### trait rate matrix
sigmasq <- 0.3145740
correlation <- (-0.632620487603683)
correlation.matrix <- diag(2)
correlation.matrix[1,2] <- correlation
correlation.matrix[2,1] <- correlation
trait.rate.matrix <- sigmasq * correlation.matrix # 0.314574, -0.199005957267441, -0.199005957267441, 0.314574

#### determinant
det.trait.rate.matrix <- log(det(trait.rate.matrix)) # 2.824245359292378
det.inv.traitRate.matrix <- det(solve(trait.rate.matrix)) # 16.84822583043347

### (3) testNonShrinkageOperationsMultipleRates
dat <- matrix(c(-2.62762948691895, -0.764018322006132,
                -1.50846427625826, -1.02686498716963,
                -0.226074849617958, -1.73165056392106), nrow = 3, ncol = 2, byrow = TRUE)

#### trait rate matrix
sigmasq <- c(0.3, 0.2)
correlation <- (-0.720107524122507)
trait.rate.matrix <- diag(2)
trait.rate.matrix[1,1] <- sigmasq[1]
trait.rate.matrix[2,2] <- sigmasq[2]
trait.rate.matrix[1,2] <- sqrt(sigmasq[1]) * sqrt(sigmasq[2]) * correlation
trait.rate.matrix[2,1] <- sqrt(sigmasq[1]) * sqrt(sigmasq[2]) * correlation
trait.rate.matrix <- sigmasq * correlation.matrix # 0.3, -0.176389599403907, -0.176389599403907, 0.2

#### determinant
det.trait.rate.matrix <- log(det(trait.rate.matrix)) # -3.544373678152544
det.inv.traitRate.matrix <- det(solve(trait.rate.matrix)) # 34.61799654333531


