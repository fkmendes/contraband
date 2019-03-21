library(ggplot2)
library(gtools)
library(gridExtra)

source("calibrated_validation_utils.R")
load(rdata.path)

n.sim <- 100
n.param <- 2
sigma.rate <- 5

param.names <- c("sigmasq", "mu")
beast.param.names <- c("BMSigmaSq", "BMMean")
mle.param.names <- c("sigmasq.mle", "mu.mle")

## template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_template.xml"
## res.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_results/"
## res.files <- mixedsort(paste0(res.path,list.files(res.path)))
## rdata.path <- gsub("_template.xml", ".RData", template.path)

## for non-ultrametric analysis (comment out or in w.r.t. to lines above)
template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_nonultra_template.xml"
res.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_nonultra_results/"
res.files <- mixedsort(paste0(res.path,list.files(res.path)))

############# DOING STUFF #############

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
table((full.df$mu >= full.df$lower.mu) & (full.df$mu <= full.df$upper.mu))
table((full.df$sigmasq >= full.df$lower.sigmasq) & (full.df$sigmasq <= full.df$upper.sigmasq))

### PLOTS ###
min.x <- min(full.df$sigmasq)
max.x <- max(full.df$sigmasq)
min.y <- min(full.df$mean.sigmasq)
max.y <- max(full.df$mean.sigmasq)
x.lab <- expression(sigma^2)
prior.mean <- 1/sigma.rate

pl1 <- get.plot("sigmasq", "mean.sigmasq", min.x, max.x, min.y, max.y, x.lab, prior.mean, full.df)
pl1

x.lab <- expression(mu)
prior.mean <- x0.mean
min.x <- min(full.df$mu)
max.x <- max(full.df$mu)
min.y <- min(full.df$mean.mu)
max.y <- max(full.df$mean.mu)

pl2 <- get.plot("mu", "mean.mu", min.x, max.x, min.y, max.y, x.lab, prior.mean, full.df)
pl2

sub.df <- full.df[full.df$"sigmasq"<.1,]
min.x <- min(sub.df$mu)
max.x <- max(sub.df$mu)
min.y <- min(sub.df$mean.mu)
max.y <- max(sub.df$mean.mu)

pl2 <- get.plot("mu", "mean.mu", min.x, max.x, min.y, max.y, x.lab, prior.mean, sub.df)
pl2

grid.arrange(pl1, pl2, ncol=2)
