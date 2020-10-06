# This R script gives us the expected values for JUnit tests
# 'PruneLikelihoodProcessTest'
# (1) 'testOneRateOnlyWithoutShrinkage'
# (2) 'testMultipleRatesWithoutShrinkage'
# (3) 'testOneRateOnlyWithShrinkage'

library(mvMORPH)
library(PCMBase)

## 'PruneLikelihoodProcessTest'

### (1) testOneRateOnlyWithoutShrinkage
tr <- read.tree(text = "((A:23.0058179,B:23.0058179):14.350951,C:37.3567689);")

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
inverse.trait.rate.matrix <- solve(trait.rate.matrix)
det.inverse.trait.rate.matrix <- det(inverse.trait.rate.matrix)

rates <- rep(1, 4)
times <- c(23.0058179, 23.0058179, 37.3567689, 14.350951)

vA = rates[1] * times[1]
vB = rates[2] * times[2]
vC = rates[3] * times[3] 
vD = rates[4] * times[4]

mA = dat[1,]
mB = dat[2,]
mC = dat[3,]

aMat.A = -0.5 * (1 / vA) # -0.0217336328651024
cMat.A = -0.5 * (1 / vA) # -0.0217336328651024
eMat.A = (1 / vA) # 0.0434672657302047
f.A = - 0.5 * log(det.trait.rate.matrix * (vA^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.56150152287921
lMat.A = cMat.A # -0.0217336328651024
mVec.A =  eMat.A * inverse.trait.rate.matrix %*% mA # -0.833127880782627, -0.743015449643922
r.A = aMat.A * t(mA) %*% inverse.trait.rate.matrix %*% mA + f.A # -5.2367146815829

aMat.B = -0.5 * (1 / vB) # -0.0217336328651024
cMat.B = -0.5 * (1 / vB) # -0.0217336328651024
eMat.B = (1 / vB) # 0.0434672657302047
f.B = - 0.5 * log(det.trait.rate.matrix * (vB^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.56150152287921
lMat.B = cMat.B # -0.0217336328651024
mVec.B = eMat.B * inverse.trait.rate.matrix %*% mB # -0.57994793022761, -0.58725740810726
r.B = aMat.B * t(mB) %*% inverse.trait.rate.matrix %*% mB + f.B # -4.46720421241219

lMat.AB = lMat.A + lMat.B
mVec.AB = mVec.A + mVec.B
r.AB = r.A + r.B

aMat.C = -0.5 * (1 / vC) # -0.0133844552064566
cMat.C = -0.5 * (1 / vC) # -0.0133844552064566
eMat.C = (1 / vC) # 0.0267689104129131
f.C = - 0.5 * log(det.trait.rate.matrix  * (vC^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -4.04626851083767
lMat.C = cMat.C # -0.0133844552064566
mVec.C = eMat.C * inverse.trait.rate.matrix %*% mC # -0.22145452354273, -0.319649013357964
r.C = aMat.C * t(mC) %*% inverse.trait.rate.matrix %*% mC + f.C # -4.40853145593445

aMat.D = -0.5 * (1 / vD) # -0.0348408966067824
cMat.D = -0.5 * (1 / vD) # -0.0348408966067824
eMat.D = (1 / vD) # 0.0696817932135647
f.D = - 0.5 * log(det.trait.rate.matrix  * (vD^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.08957059854991
lMat.D =  cMat.D - 0.25 * eMat.D * (1/(aMat.D + lMat.AB)) * t(eMat.D) # -0.0193394719769881
mVec.D = - 0.5 * eMat.D * (1/(aMat.D + lMat.AB)) * mVec.AB # -0.628706213499009, -0.591865492849479
r.D = f.D + r.AB + ((trait.nr / 2) * log(2 * pi)) - 0.5 * trait.nr * log(- 2 * (aMat.D + lMat.AB)) - 0.5 * log(det.inverse.trait.rate.matrix) - 0.25 * (1/(aMat.D + lMat.AB)) * (t(mVec.AB) %*% trait.rate.matrix %*%  mVec.AB) # -9.11979587278366

mE = c(-1.31465955080609, -1.79611255566212)
#vD_ = vD + (vA * vB / (vA + vB))
#mD_ = (vB * mA + vA * mB) / (vA + vB) 
#mE = (vC * mD_ + vD_ * mC) / (vC + vD_)

lMat.E = lMat.C + lMat.D # -0.0327239271834447
mVec.E =mVec.C + mVec.D # -0.850160737041738, -0.911514506207443
r.E = r.C + r.D # -13.5283273287181
# L = lMat.E * t(mE) %*% inverse.trait.rate.matrix %*% mE + t(mE) %*% mVec.E + r.E
# -> -12.1509000377483

# modelBM <- PCM("BM", k = trait.nr, regimes = 1)
# modelBM$X0[] = mE
# modelBM$Sigma_x[,,1] = t(chol(trait.rate.matrix))
# modelBM$Sigmae_x <- NULL
# PCMLik(t(dat),tr, modelBM)
# -> -12.1509000377483

# mvLL(tr, data = dat, method = "pic", param = list(estim = FALSE, sigma = trait.rate.matrix))
# -> -12.1509000377483

# vCD = (vC * vD_) / (vC+ vD_)
# root.2.subtract = (- 0.5 * log(det.trait.rate.matrix * (vCD^trait.nr)) - (trait.nr / 2) * log(2 * pi)) 
# L - root.2.subtract = -8.99864408103911

### (2) testMultipleRatesWithoutShrinkage
tr <- read.tree(text = "((A:23.0058179,B:23.0058179):14.350951,C:37.3567689);")

trait.nr <- 2
dat <- matrix(c(-2.62762948691895, -0.764018322006132,
                -1.50846427625826, -1.02686498716963,
                -0.226074849617958, -1.73165056392106), nrow = 3, ncol = trait.nr, byrow = TRUE)

sigmasq <- c(0.3, 0.2)
correlation <- (-0.720107524122507)
trait.rate.matrix <- diag(sigmasq)
trait.rate.matrix[1,2] <- sqrt(sigmasq[1]) * sqrt(sigmasq[2]) * correlation
trait.rate.matrix[2,1] <- sqrt(sigmasq[1]) * sqrt(sigmasq[2]) * correlation

det.trait.rate.matrix <- det(trait.rate.matrix)
inverse.trait.rate.matrix <- solve(trait.rate.matrix)
det.inverse.trait.rate.matrix <- det(inverse.trait.rate.matrix)

rates <- rep(1, 4)
times <- c(23.0058179, 23.0058179, 37.3567689, 14.350951)

vA = rates[1] * times[1]
vB = rates[2] * times[2]
vC = rates[3] * times[3] 
vD = rates[4] * times[4]

mA = dat[1,]
mB = dat[2,]
mC = dat[3,]

aMat.A = -0.5 * (1 / vA) # -0.0217336328651024
cMat.A = -0.5 * (1 / vA) # -0.0217336328651024
eMat.A = (1 / vA) # 0.0434672657302047
f.A = - 0.5 * log(det.trait.rate.matrix * (vA^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.20143736344913
lMat.A = cMat.A # -0.0217336328651024
mVec.A =  eMat.A * inverse.trait.rate.matrix %*% mA # -0.993572327994746, -1.04232806169593
r.A = aMat.A * t(mA) %*% inverse.trait.rate.matrix %*% mA + f.A # -4.90498620500039

aMat.B = -0.5 * (1 / vB) # -0.0217336328651024
cMat.B = -0.5 * (1 / vB) # -0.0217336328651024
eMat.B = (1 / vB) # 0.0434672657302047
f.B = - 0.5 * log(det.trait.rate.matrix * (vB^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.20143736344913
lMat.B = cMat.B # -0.0217336328651024
mVec.B = eMat.B * inverse.trait.rate.matrix %*% mB # -0.726524972304187, -0.863932310440079
r.B = aMat.B * t(mB) %*% inverse.trait.rate.matrix %*% mB + f.B # -4.19297676715206

lMat.AB = lMat.A + lMat.B
mVec.AB = mVec.A + mVec.B
r.AB = r.A + r.B

aMat.C = -0.5 * (1 / vC) # -0.0133844552064566
cMat.C = -0.5 * (1 / vC) # -0.0133844552064566
eMat.C = (1 / vC) # 0.0267689104129131
f.C = - 0.5 * log(det.trait.rate.matrix  * (vC^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.68620435140759
lMat.C = cMat.C # -0.0133844552064566
mVec.C = eMat.C * inverse.trait.rate.matrix %*% mC # -0.32495184010392, -0.518362618567832
r.C = aMat.C * t(mC) %*% inverse.trait.rate.matrix %*% mC + f.C # -4.17174753097916

aMat.D = -0.5 * (1 / vD) # -0.0348408966067824
cMat.D = -0.5 * (1 / vD) # -0.0348408966067824
eMat.D = (1 / vD) # 0.0696817932135647
f.D = - 0.5 * log(det.trait.rate.matrix  * (vD^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -2.72950643911983
lMat.D =  cMat.D - 0.25 * eMat.D * (1/(aMat.D + lMat.AB)) * t(eMat.D) # -0.0193394719769881
mVec.D = - 0.5 * eMat.D * (1/(aMat.D + lMat.AB)) * mVec.AB # -0.765306328290814, -0.84813407120176
r.D = f.D + r.AB + ((trait.nr / 2) * log(2 * pi)) - 0.5 * trait.nr * log(- 2 * (aMat.D + lMat.AB)) - 0.5 * log(det.inverse.trait.rate.matrix) - 0.25 * (1/(aMat.D + lMat.AB)) * (t(mVec.AB) %*% trait.rate.matrix %*%  mVec.AB) # -8.44680149853574

mE = c(-1.31465955080609, -1.2374605274288)
#vD_ = vD + (vA * vB / (vA + vB))
#mD_ = (vB * mA + vA * mB) / (vA + vB) 
#mE = (vC * mD_ + vD_ * mC) / (vC + vD_)

lMat.E = lMat.C + lMat.D # -0.0327239271834447
mVec.E =mVec.C + mVec.D # -1.09025816839473, -1.36649668976959
r.E = r.C + r.D # -12.6185490295149
# L = lMat.E * t(mE) %*% inverse.trait.rate.matrix %*% mE + t(mE) %*% mVec.E + r.E
# -> -11.0563970153267

# modelBM <- PCM("BM", k = trait.nr, regimes = 1)
# modelBM$X0[] = mE
# modelBM$Sigma_x[,,1] = t(chol(trait.rate.matrix))
# modelBM$Sigmae_x <- NULL
# PCMLik(t(dat),tr, modelBM)
# -> -11.0563970153267

# mvLL(tr, data = dat, method = "pic", param = list(estim = FALSE, sigma = trait.rate.matrix))
# -> -11.0563970153267

# vCD = (vC * vD_) / (vC+ vD_)
# root.2.subtract = (- 0.5 * log(det.trait.rate.matrix * (vCD^trait.nr)) - (trait.nr / 2) * log(2 * pi)) 
# L - root.2.subtract = -8.26420521804756

### (3) testOneRateOnlyWithShrinkage
tr <- read.tree(text = "((A:12.4420263,B:12.4420263):42.9258211,C:43.5702874);")

trait.nr <- 2
dat <- matrix(c(1.0, 2.0,
                3.0, 5.0,
                2.0, 4.0), nrow = 3, ncol = trait.nr, byrow = TRUE)

lambda <- estimate.lambda(dat) # 0.25925925925926

shrinkage.correlation.matrix <- matrix(cor.shrink(dat, verbose = FALSE), trait.nr, trait.nr)
sigmasq <- 0.1543038
trait.rate.matrix <- sigmasq * shrinkage.correlation.matrix

lMat = chol(trait.rate.matrix) 
# inverse of lMat -> A = L.inverse
aMat = solve(lMat)
# transformed data -> Z = M * A
transformed.trait.values = dat %*% aMat

det.trait.rate.matrix <- det(trait.rate.matrix) 
inverse.trait.rate.matrix <- solve(trait.rate.matrix)
det.inverse.trait.rate.matrix <- det(inverse.trait.rate.matrix) 

rates <- rep(1, 4)
times = c(12.4420263, 12.4420263, 43.5702874, 42.9258211)

vA = rates[1] * times[1]
vB = rates[2] * times[2]
vC = rates[3] * times[3] 
vD = rates[4] * times[4]

mA = transformed.trait.values[1,]
mB = transformed.trait.values[2,]
mC = transformed.trait.values[3,]

aMat.A = -0.5 * (1 / vA) # -0.0401863802522263
cMat.A = -0.5 * (1 / vA) # -0.0401863802522263
eMat.A = (1 / vA) # 0.0803727605044526
f.A = - 0.5 * log(det.trait.rate.matrix * (vA^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -2.11356981107731
lMat.A = cMat.A # -0.0401863802522263
mVec.A = eMat.A %*% mA # 0.204607040785534, 0.379446710223391
r.A = aMat.A * t(mA) %*% mA + f.A # -3.26970682734958

aMat.B = -0.5 * (1 / vB) # -0.0401863802522263
cMat.B = -0.5 * (1 / vB) # -0.0401863802522263
eMat.B = (1 / vB) # 0.0803727605044526
f.B = - 0.5 * log(det.trait.rate.matrix * (vB^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -2.11356981107731
lMat.B = cMat.B # -0.0401863802522263
mVec.B = eMat.B %*% mB # 0.613821122356603, 0.840175260824958
r.B = aMat.B * t(mB) %*% mB + f.B # -8.84887933857218

lMat.AB = lMat.A + lMat.B
mVec.AB = mVec.A + mVec.B
r.AB = r.A + r.B

aMat.C = -0.5 * (1 / vC) # -0.0114757103943271
cMat.C = -0.5 * (1 / vC) # -0.0114757103943271
eMat.C = (1 / vC)  # 0.0229514207886542
f.C = - 0.5 * log(det.trait.rate.matrix * (vC^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.36686528756809
lMat.C = cMat.C # -0.0114757103943271
mVec.C = eMat.C %*% mC # 0.116856065659957, 0.216711260346101
r.C = aMat.C * t(mC) %*% mC + f.C # -4.68746131952019

aMat.D = -0.5 * (1 / vD) # -0.0116480008346305
cMat.D = -0.5 * (1 / vD) # -0.0116480008346305
eMat.D = (1 / vD) # 0.023296001669261
f.D = - 0.5 * log(det.trait.rate.matrix * (vD^trait.nr)) - (trait.nr / 2) * log(2 * pi) # -3.3519633864921
lMat.D = cMat.D - 0.25 * eMat.D %*% solve(aMat.D + lMat.AB) %*% t(eMat.D) # -0.0101735952606144
mVec.D = - 0.5 * eMat.D %*% solve(aMat.D + lMat.AB) %*% mVec.AB # 0.10359675130525, 0.154379919596158
r.D = f.D + r.AB + ((trait.nr / 2) * log(2 * pi)) - 0.5 * log(((-2 * (aMat.D + lMat.AB))^trait.nr) * det.inverse.trait.rate.matrix) - 0.25 * solve(aMat.D + lMat.AB) * (mVec.AB %*%  t(mVec.AB)) # -8.32455362290912

lMat.E = lMat.C + lMat.D # -0.0216493056549415
mVec.E = mVec.C + mVec.D # 0.220452816965207, 0.371091179942258
r.E = r.C + r.D # -13.0120149424293
# L = lMat.E * t(mE) %*% mE + t(mE) %*% t(mVec.E) + r.E 
# -> -10.86058210541

vD_ = vD + (vA * vB / (vA + vB))
mD_ = (vB * mA + vA * mB) / (vA + vB) 
mE = (vC * mD_ + vD_ * mC) / (vC + vD_) # -> 5.09145236523758, 8.57050997054855

vCD = (vC * vD_) / (vC+ vD_)
root.2.subtract = (- 0.5 * log(det.trait.rate.matrix * (vCD^trait.nr)) - (trait.nr / 2) * log(2 * pi)) 
L = lMat.E * t(mE) %*% mE + t(mE) %*% t(mVec.E) + r.E - root.2.subtract
# -> -8.12845753837811

# res1 = mvLL(tr, data = dat, method = "pic", param = list(estim = FALSE, sigma = trait.rate.matrix))
# res1$logl 
# -> -10.86058210541
# res1$theta
# -> 2, 3.76503645376053

# modelBM <- PCM("BM", k = trait.nr, regimes = 1)
# modelBM$X0[] = c(2, 3.76503645376053)
# modelBM$Sigma_x[,,1] = t(chol(trait.rate.matrix))
# modelBM$Sigmae_x <- NULL
# res2 = PCMLik(t(dat),tr, modelBM)
# res2[1]
# -> -10.86058210541



