library(ggplot2)
library(gtools)
library(sjPlot)
library(ape)
library(phytools)
source("calibrated_validation_utils.R")

args = commandArgs(trailingOnly=TRUE)

### SCRIPT FLAGS AND PATH VARIABLES ###

cal.validation.folder <- args[1]
rdata.path <- args[2]
n.sim <- as.numeric(args[3])
n.sim <- 100
n.param <- 6
sigma.rate <- 5
rv.mean <- 0.0 # root value
th.mean <- 0.0 # theta
alpha.mean <- 1.0
job.prefix <- args[4] # e.g., "OUThreeOptMVN"
tree.type <- args[5] # e.g., "ultra" or "nonultra"

param.labs <- c(expression(sigma^2), expression(y[0]), expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(alpha))
param.names <- c("sigmasq", "rv", "theta1", "theta2", "theta3", "alpha")
beast.param.names <- c("OUSigmaSq", "OURootValue", "OUTheta1", "OUTheta3", "OUTheta2", "OUAlpha")
mle.param.names <- c("sigmasq.mle", "rv.mle", "theta1.mle", "theta2.mle", "theta3.mle", "alpha.mle")
prior.means <- c(1/sigma.rate, rv.mean, th.mean, th.mean, th.mean, alpha.mean)

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
table((full.df$sigmasq >= full.df$lower.sigmasq) & (full.df$sigmasq <= full.df$upper.sigmasq))
table((full.df$rv >= full.df$lower.rv) & (full.df$rv <= full.df$upper.rv))
table((full.df$theta1 >= full.df$lower.theta1) & (full.df$theta1 <= full.df$upper.theta1))
table((full.df$theta2 >= full.df$lower.theta2) & (full.df$theta2 <= full.df$upper.theta2))
table((full.df$theta3 >= full.df$lower.theta3) & (full.df$theta3 <= full.df$upper.theta3))
table((full.df$alpha >= full.df$lower.alpha) & (full.df$alpha <= full.df$upper.alpha))

## ultrametric
## FALSE  TRUE
##     5    95

## FALSE  TRUE
##     3    97

## FALSE  TRUE
##     6    94

## FALSE  TRUE
##     2    98

## FALSE  TRUE
##    13    87

## FALSE  TRUE
##     8    92

## nonultrametric

all.plots <- vector("list", n.param)
for (i in 1:n.param) {
    x.lab = param.labs[i]
    min.x = min(full.df[,param.names[i]])
    max.x = max(full.df[,param.names[i]])
    min.y = min(full.df[,paste0("mean.",param.names[i])])
    max.y = max(full.df[,paste0("mean.",param.names[i])])
    print(param.names[i]);     print(paste0("mean.",param.names[i]))
    all.plots[[i]] = get.plot(param.names[i], paste0("mean.",param.names[i]),
                              min.x, max.x, min.x, max.x, x.lab, prior.means[i],
                              full.df)
    names(all.plots)[i] = paste0("plot", i)
}
list2env(all.plots, .GlobalEnv) # sending plots in list into environment so I cna use plot_grid

png(paste0(cal.validation.folder, job.prefix, "_", tree.type, "_graphs.png"), height=30, width=20, unit="cm", res=300)
plot_grid(mget(paste0("plot", 1:n.param)))
dev.off()

png(paste0(cal.validation.folder, job.prefix, "_", tree.type, "_tree.png"), height=10, width=10, unit="cm", res=300)
plotSimmap(tr, ftype="off")
dev.off()
