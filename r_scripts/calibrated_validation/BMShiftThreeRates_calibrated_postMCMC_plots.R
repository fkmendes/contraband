library(ggplot2)
library(gtools)
library(sjPlot)
library(ape)
library(gridExtra)
library(phytools)
source("calibrated_validation_utils.R")

args = commandArgs(trailingOnly=TRUE)

### SCRIPT FLAGS AND PATH VARIABLES ###

cal.validation.folder <- args[1]
rdata.path <- args[2]
n.sim <- 100
n.sim <- as.numeric(args[3])
n.param <- 4
sigma.rate <- 10
x0s.mean <- 0
job.prefix <- args[4] # e.g., "BMMVNShiftThreeRates" or "BMPruneShiftThreeRates"
tree.type <- args[5] # e.g., "ultrametric" or "nonultrametric"

param.labs <- c(expression(sigma[1]^2), expression(sigma[2]^2), expression(sigma[3]^2), expression(mu))
param.names <- c("sigma1", "sigma2", "sigma3", "mu")
beast.param.names <- c("rateValues1", "rateValues3", "rateValues2", "BMMean")
mle.param.names <- c("sigma1.mle", "sigma2.mle", "sigma3.mle", "mu.mle")
prior.means <- c(1/sigma.rate, 1/sigma.rate, 1/sigma.rate, x0s.mean)

res.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_results/")
res.files <- mixedsort(paste0(res.path,list.files(res.path)))
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
table((full.df$mu >= full.df$lower.mu) & (full.df$mu <= full.df$upper.mu))
table((full.df$sigma1 >= full.df$lower.sigma1) & (full.df$sigma1 <= full.df$upper.sigma1))
table((full.df$sigma2 >= full.df$lower.sigma2) & (full.df$sigma2 <= full.df$upper.sigma2))
table((full.df$sigma3 >= full.df$lower.sigma3) & (full.df$sigma3 <= full.df$upper.sigma3))

## FALSE  TRUE
##     7    93

## FALSE  TRUE
##     4    96

## FALSE  TRUE
##     1    99

## FALSE  TRUE
##     5    95

## FALSE  TRUE
##    11    89

## FALSE  TRUE
##     3    97

## FALSE  TRUE
##     3    97

## FALSE  TRUE
##     4    96

### PLOTS ###

all.plots <- vector("list", n.param)
for (i in 1:n.param) {
    x.lab = param.labs[i]
    min.x = min(full.df[,param.names[i]])
    max.x = max(full.df[,param.names[i]])
    min.y = min(full.df[,paste0("mean.",param.names[i])])
    max.y = max(full.df[,paste0("mean.",param.names[i])])
    all.plots[[i]] = get.plot(param.names[i], paste0("mean.",param.names[i]),
                              min.x, max.x, min.y, max.y, x.lab, prior.means[i],
                              full.df)
    names(all.plots)[i] = paste0("plot", i)
}
list2env(all.plots, .GlobalEnv) # sending plots in list into environment so I cna use plot_grid

png(paste0(cal.validation.folder, job.prefix, "_", tree.type, "_graphs.png"), height=20, width=20, unit="cm", res=300)
plot_grid(mget(paste0("plot", 1:4)))
dev.off()

png(paste0(cal.validation.folder, job.prefix, "_", tree.type, "_tree.png"), height=10, width=10, unit="cm", res=300)
plotSimmap(tr, ftype="off")
dev.off()
