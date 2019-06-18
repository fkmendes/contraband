library(ggplot2)
library(gtools)
library(ggpubr)
library(ape)
source("calibrated_validation_utils.R")

args = commandArgs(trailingOnly=TRUE)

### SCRIPT FLAGS AND PATH VARIABLES ###

cal.validation.folder <- args[1]
## cal.validation.folder <- "./"
rdata.path <- args[2]
## rdata.path <- "BMMVNShiftLikelihoodOneRateFBDOneTrait_nonultra.RData"
n.sim <- as.numeric(args[3])
## n.sim <- 100
job.prefix <- args[4]
## job.prefix <- "BMMVNShiftOneRateFBD"
n.param <- as.numeric(args[5])
## n.param <- 2

param.names <- strsplit(args[6], ",")[[1]]
## param.names <- c("sigmasq","mu")
beast.param.names <- strsplit(args[7], ",")[[1]]
## beast.param.names <- c("rateValues", "BMMean")
param.labs.preparse <- strsplit(args[8], ",")[[1]]
prior.means.preparse <- strsplit(args[9], ",")[[1]]
param.labs <- prior.means <- c()
for (i in 1:n.param) {
    param.labs[i] = eval(parse(text=param.labs.preparse[i]))
    prior.means[i] = eval(parse(text=prior.means.preparse[i]))
}
mle.param.names <- strsplit(args[11], ",")[[1]]

## param.labs <- c(expression(sigma^2), expression(mu))
## param.names <- c("sigmasq", "mu")
## beast.param.names <- c("rateValues", "BMMean")
## mle.param.names <- c("sigmasq.mle", "mu.mle")
## prior.means <- c(0.2, 0.0)

tree.type <- "nonultra"

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
table((full.df$sigmasq >= full.df$lower.sigmasq) & (full.df$sigmasq <= full.df$upper.sigmasq))

## nonultrametric
## mu
##   FALSE  TRUE
##       3    97
## sigma
##   FALSE  TRUE
##      8    92

### PLOTS ###

all.plots <- vector("list", n.param)
for (i in 1:n.param) {
    x.lab = param.labs[i]
    min.x = min(full.df[,param.names[i]])
    max.x = max(full.df[,param.names[i]])
    min.y = min(full.df[,paste0("mean.",param.names[i])], na.rm=TRUE)
    max.y = max(full.df[,paste0("mean.",param.names[i])], na.rm=TRUE)
    all.plots[[i]] = get.plot(param.names[i], paste0("mean.",param.names[i]),
                              min.x, max.x, min.y, max.y, x.lab, prior.means[i],
                              full.df)
    names(all.plots)[i] = paste0("plot", i)
}
list2env(all.plots, .GlobalEnv) # sending plots in list into environment so I cna use plot_grid

png(paste0(cal.validation.folder, job.prefix, "_", tree.type, "_graphs.png"), height=16, width=10, unit="cm", res=300)
ggarrange(plotlist=all.plots, ncol=1, nrow=2)
dev.off()
