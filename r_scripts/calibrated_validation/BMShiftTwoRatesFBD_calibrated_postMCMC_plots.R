library(ggplot2)
library(gtools)
library(ggpubr)
library(ape)
source("calibrated_validation_utils.R")

args = commandArgs(trailingOnly=TRUE)

### SCRIPT FLAGS AND PATH VARIABLES ###

tree.type <- "nonultra"

cal.validation.folder <- "./"
rdata.path <- "BMMVNShiftLikelihoodTwoRatesFBDfixedOneTrait_nonultra.RData"
n.sim <- 100
job.prefix <- "BMMVNShiftTwoRatesFBDfixed"
n.param <- 3
param.names <- c("sigma1", "sigma2", "mu")
beast.param.names <- c("rateValues1", "rateValues2", "BMMean")
param.labs <- c(expression(sigma[1]^2), expression(sigma[2]^2), expression(y[0]))
mle.param.names <- c("sigma1.mle", "sigma2.mle", "mu.mle")
prior.means <- c(1.516004, 1.516004, 0.0)

## cal.validation.folder <- args[1]
## rdata.path <- args[2]
## n.sim <- as.numeric(args[3])
## job.prefix <- args[4]
## n.param <- as.numeric(args[5])
## param.names <- strsplit(args[6], ",")[[1]]
## beast.param.names <- strsplit(args[7], ",")[[1]]
## param.labs.preparse <- strsplit(args[8], ",")[[1]]
## prior.means.preparse <- strsplit(args[9], ",")[[1]]
## param.labs <- prior.means <- c()
## for (i in 1:n.param) {
##     param.labs[i] = eval(parse(text=param.labs.preparse[i]))
##     prior.means[i] = eval(parse(text=prior.means.preparse[i]))
## }
## mle.param.names <- strsplit(args[11], ",")[[1]]

res.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_results/")
res.files <- mixedsort(paste0(res.path,list.files(res.path, pattern = "[0-9]\\.log$")))
rate.files <- mixedsort(paste0(res.path,list.files(res.path, pattern = "\\.trees\\.log$")))
load(rdata.path)

############# DOING STUFF #############

n.cols = n.param * 3
log.df <- data.frame(matrix(ncol=n.cols, nrow=n.sim))
names(log.df) <- as.vector(outer(c("lower", "upper", "mean"), param.names, paste, sep="."))
cols <- seq(length.out=n.param, by=3) # 3=lower, upper, mean
cols.2.compare <- c(3, 6)
for (i in 1:n.sim) {
    this.sim.rate.df = read.table(rate.files[i], header=TRUE, row.names=1, sep="\t")
    names(this.sim.rate.df) = beast.param.names[1:2]
    this.sim.df = read.delim(res.files[i], comment.char='#')
    this.sim.df = cbind.fill(this.sim.rate.df, this.sim.df)

    k = 1
    for (j in cols) {
        beast.param.name = as.character(beast.param.names[k])
        ## log.df[i,j:(j+2)] = get.95(this.sim.df[102:1001,beast.param.name])
        log.df[i,j:(j+2)] = get.95(this.sim.df[102:nrow(this.sim.df),beast.param.name]) # for incomplete runs
        k = k+1
    }
}

# putting true values and estimated ones together
full.df <- cbind(true.param.df, log.df)
bool.vec <- full.df[,1] > full.df[,2]
cols <- c((2*n.param+1):(2*n.param+3))
for (i in cols) {
    full.df <- flip.w.bool(full.df, bool.vec, i, i+3)
}
## cols.2.flip <- c(1, 4) # 1 and 2 should flip, 4 and 5 should flip
## full.df <- flip.trace.var(full.df, cols.2.flip)

plot.hdi <- vector("list", n.param)
for (i in 1:n.param) {
    upper = paste0("upper.", param.names[i])
    lower = paste0("lower.", param.names[i])
    plot.hdi[[i]] = (full.df[,upper] < full.df[,param.names[i]] | full.df[,lower] > full.df[,param.names[i]])
}

# coverage
table((full.df$mu >= full.df$lower.mu) & (full.df$mu <= full.df$upper.mu))
table((full.df$sigma1 >= full.df$lower.sigma1) & (full.df$sigma1 <= full.df$upper.sigma1))
table((full.df$sigma2 >= full.df$lower.sigma2) & (full.df$sigma2 <= full.df$upper.sigma2))

## nonultrametric
## mu
##   FALSE  TRUE
##       6    91
## sigma1
##   FALSE  TRUE
##      6    92
## sigma2
##   FALSE  TRUE
##      12    90

### PLOTS ###

all.plots <- vector("list", n.param)
for (i in 1:n.param) {
    x.lab = param.labs[i]
    min.x = min(full.df[,param.names[i]])
    max.x = max(full.df[,param.names[i]])
    min.y = min(full.df[,paste0("mean.",param.names[i])])
    max.y = max(full.df[,paste0("mean.",param.names[i])])
    all.plots[[i]] = get.plot(param.names[i], paste0("mean.",param.names[i]),
                              min.x, max.x, min.x, max.x, x.lab, prior.means[i],
                              full.df, plot.hdi[[i]])
    names(all.plots)[i] = paste0("plot", i)
}
list2env(all.plots, .GlobalEnv) # sending plots in list into environment so I cna use plot_grid

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_sigsq1.pdf"), height=2, width=2.5)
plot1
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_sigsq2.pdf"), height=2, width=2.5)
plot2
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_y0.pdf"), height=2, width=2.5)
plot3
dev.off()

## png(paste0(paste0(cal.validation.folder, "figs/"), job.prefix, "_", tree.type, "_graphs.png"), height=24, width=10, unit="cm", res=300)
## ggarrange(plotlist=all.plots, ncol=1, nrow=3)
## dev.off()
