library(ggplot2)
library(gtools)
library(gridExtra)
library(sjPlot)

source("calibrated_validation_utils.R")

n.sim <- 100
n.param <- 4
sigma.rate <- 5
rv.mean <- 0.0 # root value
th.mean <- 0.0 # theta
alpha.mean <- 1.0

param.labs <- c(expression(sigma^2), expression(mu), expression(theta), expression(alpha))
param.names <- c("sigmasq", "rv", "theta", "alpha")
beast.param.names <- c("OUSigmaSq", "OURootValue", "OUTheta", "OUAlpha")
mle.param.names <- c("sigmasq.mle", "rv.mle", "theta.mle", "alpha.mle")
prior.means <- c(1/sigma.rate, rv.mean, th.mean, alpha.mean)

template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/OUVanillaMVNLikelihoodOneTrait_fixedtree_template.xml"
res.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/OUVanillaMVNOneTrait_results/"
res.files <- mixedsort(paste0(res.path,list.files(res.path)))
rdata.path <- gsub("_template.xml", ".RData", template.path)
load(rdata.path)

############# DOING STUFF #############

n.cols = n.param * 3
log.df <- data.frame(matrix(ncol=n.cols, nrow=n.sim))
names(log.df) <- as.vector(outer(c("lower", "upper", "mean"), param.names, paste, sep="."))
cols <- seq(length.out=n.param, by=3) # 3=lower, upper, mean
for (i in 1:n.sim) {
    this.sim.df <- read.delim(res.files[i], comment.char='#')
    
    k = 1
    for (j in cols) {
        beast.param.name = as.character(beast.param.names[k])
        log.df[i,j:(j+2)] = get.95(this.sim.df[102:1001,beast.param.name])
        k = k+1
    }
}

# putting true values and estimated ones together
full.df <- cbind(true.param.df, log.df)

# coverage
table((full.df$sigmasq >= full.df$lower.sigmasq) & (full.df$sigmasq <= full.df$upper.sigmasq))
table((full.df$rv >= full.df$lower.rv) & (full.df$rv <= full.df$upper.rv))
