# This R script gives us the expected values for JUnit tests
# 'BMPruneLikelihoodUtils'
# (1) 'testPopulateTraitValues'
# (2) 'testPopulatePCMParams'
# (3) 'testPopulatePCMParamsWithShrinkage'
# (4) 'testPopulateEstimatedVarianceMatrix'
# (5) 'testPopulateGivenPopulationVariance'

## 'BMPruneLikelihoodUtils'

### (1) testPopulateTraitValues
dat <- matrix(c(0.983714690867666, -7.54729477473779,
                -7.86424514338822, -2.97908131550921,
                7.23079460908758, -0.780647498381348,
                -1.39605330265115, 3.72028693114977), nrow = 4, ncol = 2, byrow = TRUE)

trait.values.array <- c()
species.nr <- c(3, 4, 1, 2)
for(i in species.nr ) {
  trait.values.array  = c(trait.values.array, paste0(dat[i,], collapse = ", "))
}
(paste0(trait.values.array, collapse = ", ")) # 7.23079460908758, -0.780647498381348, -1.39605330265115, 3.72028693114977, 0.983714690867666, -7.54729477473779, -7.86424514338822, -2.97908131550921

# the first species
paste0(dat[which(species.nr == 1),], collapse = ", ") # 7.23079460908758, -0.780647498381348
# the second species
paste0(dat[which(species.nr == 2),], collapse = ", ") # -1.39605330265115, 3.72028693114977
# the third species
paste0(dat[which(species.nr == 3),], collapse = ", ") # 0.983714690867666, -7.54729477473779
# the fourth species
paste0(dat[which(species.nr == 4),], collapse = ", ") # -7.86424514338822, -2.97908131550921

### (2) testPopulatePCMParams
trait.nr <- 2
dat <- matrix(c(-2.62762948691895, -1.56292164859448,
                -1.50846427625826, -1.59482814741543,
                -0.226074849617958, -2.11000367246907), nrow = 3, ncol = trait.nr, byrow = TRUE)

sigmasq <- 0.3145740
correlation <- (-0.632620487603683)
correlation.matrix <- diag(trait.nr)
correlation.matrix[1,2] <- correlation
correlation.matrix[2,1] <- correlation
trait.rate.matrix <- sigmasq * correlation.matrix

det.trait.rate.matrix <- det(trait.rate.matrix) 
inv.trait.rate.matrix <- solve(trait.rate.matrix)
det.inv.traitRate.matrix <- det(inv.trait.rate.matrix) 

times <- c(23.0058179, 23.0058179, 37.3567689, 14.350951)
vA = times[1]
vB = times[2]
vC = times[3] 
vD = times[4]

mA = dat[1,]
mB = dat[2,]
mC = dat[3,]  

aMat.A = -0.5 * (1 / vA) # -0.02173363286510235
cMat.A = -0.5 * (1 / vA) # -0.02173363286510235
eMat.A = (1 / vA) # 0.04346726573020471
f.A = - 0.5 * log(det.trait.rate.matrix * (vA^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.561501522879213

lMat.A = cMat.A # -0.02173363286510235
mVec.A = eMat.A * inv.trait.rate.matrix %*% mA # -0.8331278807826272, -0.7430154496439217
r.A = aMat.A * t(mA) %*% inv.trait.rate.matrix %*% mA + f.A # -5.2367146815829

aMat.B = -0.5 * (1 / vB) # -0.02173363286510235
cMat.B = -0.5 * (1 / vB) # -0.02173363286510235
eMat.B = (1 / vB) # 0.04346726573020471
f.B = - 0.5 * log(det.trait.rate.matrix * (vB^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.561501522879213

lMat.B = cMat.B # -0.02173363286510235
mVec.B = eMat.B * inv.trait.rate.matrix %*% mB # -0.57994793022761, -0.58725740810726
r.B = aMat.B * t(mB) %*% inv.trait.rate.matrix %*% mB + f.B # -4.467204212412191


lMat.AB = lMat.A + lMat.B # -0.0434672657302047
mVec.AB = mVec.A + mVec.B # -1.41307581101024, -1.33027285775118
r.AB = r.A + r.B # -9.70391889399509

aMat.D = -0.5 * (1 / vD) 
cMat.D = -0.5 * (1 / vD)
eMat.D = (1 / vD) 
f.D = - 0.5 * log(det.trait.rate.matrix * (vD^trait.nr)) - (trait.nr / 2) * log(2 * pi)

lMat.D = cMat.D - 0.25 * (1 / (aMat.D + lMat.AB)) * eMat.D * eMat.D # -0.0193394719769881
mVec.D = - 0.5 * (1 / (aMat.D + lMat.AB))  *  eMat.D *  mVec.AB # -0.628706213499009, -0.591865492849479
r.D = f.D + r.AB + ((trait.nr/ 2) * log(2 * pi)) - 0.5 * log(((-2 * (aMat.D + lMat.AB))^trait.nr) * det.inv.traitRate.matrix) - 0.25 * (1 / (aMat.D + lMat.AB)) * (t(mVec.AB) %*% trait.rate.matrix %*%  mVec.AB) # -9.11979587278366
  
### (3) testPopulatePCMParamsWithShrinkage
trait.nr <- 2
dat <- matrix(c(1, 2,
                3, 5,
                2, 4), nrow = 3, ncol = trait.nr, byrow = TRUE)

shrinkage.correlation.matrix <- matrix(cor.shrink(dat, verbose = FALSE), trait.nr, trait.nr)
sigmasq <- 0.1543038
trait.rate.matrix <- sigmasq * shrinkage.correlation.matrix

lMat = chol(trait.rate.matrix) 
# inverse of lMat -> A = L.inverse
aMat = solve(lMat)
# transformed data -> Z = M * A
transformed.trait.values = dat %*% aMat

det.trait.rate.matrix <- det(trait.rate.matrix) 
inv.trait.rate.matrix <- solve(trait.rate.matrix)
det.inv.traitRate.matrix <- det(inv.trait.rate.matrix) 

times = c(12.4420263, 12.4420263, 43.5702874, 42.9258211)
vA = times[1]
vB = times[2]
vC = times[3] 
vD = times[4]

mA = transformed.trait.values[1,]
mB = transformed.trait.values[2,]
mC = transformed.trait.values[3,]  

aMat.A = -0.5 * (1 / vA) # -0.0401863802522263
cMat.A = -0.5 * (1 / vA) # -0.0401863802522263
eMat.A = (1 / vA) # 0.0803727605044526
f.A = - 0.5 * log(det.trait.rate.matrix * (vA^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.561501522879213

lMat.A = cMat.A # -0.0401863802522263
mVec.A = eMat.A *  mA # 0.0803727605044526, 0.160745521008905
r.A = aMat.A * t(mA) %*%  mA + f.A # -3.26970682734958

aMat.B = -0.5 * (1 / vB) # -0.0401863802522263
cMat.B = -0.5 * (1 / vB) # -0.0401863802522263
eMat.B = (1 / vB) # 0.0803727605044526
f.B = - 0.5 * log(det.trait.rate.matrix * (vB^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -2.11356981107731

lMat.B = cMat.B # -0.0401863802522263
mVec.B = eMat.B * mB # 0.613821122356603, 0.840175260824958
r.B = aMat.B * t(mB) %*% mB + f.B # -8.84887933857218


lMat.AB = lMat.A + lMat.B # -0.0803727605044526
mVec.AB = mVec.A + mVec.B # 0.818428163142138, 1.21962197104835
r.AB = r.A + r.B # -12.1185861659218

aMat.D = -0.5 * (1 / vD) # -0.0116480008346305
cMat.D = -0.5 * (1 / vD) # -0.0116480008346305
eMat.D = (1 / vD) # 0.023296001669261
f.D = - 0.5 * log(det.trait.rate.matrix * (vD^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.3519633864921

lMat.D = cMat.D - 0.25 * (1 / (aMat.D + lMat.AB)) * eMat.D * eMat.D # -0.0101735952606144
mVec.D = - 0.5 * (1 / (aMat.D + lMat.AB))  *  eMat.D *  mVec.AB # 0.10359675130525, 0.154379919596158
r.D = f.D + r.AB + ((trait.nr/ 2) * log(2 * pi)) - 0.5 * log(((-2 * (aMat.D + lMat.AB))^trait.nr) * det.inv.traitRate.matrix) - 0.25 * (1 / (aMat.D + lMat.AB)) * (t(mVec.AB) %*% mVec.AB) # -8.32455362290912

### (4) testPopulateEstimatedVarianceMatrix
dat <- matrix(c(-2.62762948691895, -1.56292164859448,
                -1.50846427625826, -1.59482814741543,
                -0.226074849617958, -2.11000367246907), nrow = 3, ncol = 2, byrow = TRUE)

lambda <- estimate.lambda.var(dat, verbose = FALSE)

v <- diag(var(dat))
target <- median(v)
shrinkage.variance.vector <- lambda * target + (1 - lambda) * v
# alternative: shrinkage.variance.vector <- c(var.shrink(dat, verbose = FALSE))

std <- 1 / sqrt(shrinkage.variance.vector)
transformed.dat <- dat %*% (diag(std))
# transformed.dat[1, ] -> -2.5567665115854, -2.25079276020313
# transformed.dat[2, ] -> -1.46778340122159, -2.29674190718321
# transformed.dat[3, ] -> -0.219977971587107, -3.03865583682093

### (5) testPopulateGivenPopulationVariance
dat <- matrix(c(-2.62762948691895, -1.56292164859448,
                -1.50846427625826, -1.59482814741543,
                -0.226074849617958, -2.11000367246907), nrow = 3, ncol = 2, byrow = TRUE)

population.variance <- 0.3

std <- rep(1 / sqrt(population.variance), 2)
transformed.dat <- dat %*% (diag(std))
# transformed.dat[1, ] -> -4.79737314250412, -2.8534914751612
# transformed.dat[2, ] -> -2.75406637099118, -2.91174450561202
# transformed.dat[3, ] -> -0.412754316067813, -3.85232202610017


