library(ggplot2)
library(gtools)

source("calibrated_validation_utils.R")

n.sim <- 100
n.param <- 2
param.names <- c("sigmasq", "mu")
beast.param.names <- c("BMSigmaSq", "BMMean")

template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_template.xml"
res.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_results/"
res.files <- mixedsort(paste0(res.path,list.files(res.path)))
rdata.path <- gsub("_template.xml", ".RData", template.path)
load(rdata.path)

### DOING STUFF ###

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

full.df <- cbind(true.param.df, log.df)

table((full.df$mu >= full.df$lower.mu) & (full.df$mu <= full.df$upper.mu))
table((full.df$sigmasq >= full.df$lower.sigmasq) & (full.df$sigmasq <= full.df$upper.sigmasq))
