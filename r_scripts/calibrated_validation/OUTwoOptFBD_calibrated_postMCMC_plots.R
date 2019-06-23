library(ggplot2)
library(gtools)
library(ggpubr)
library(ape)
source("calibrated_validation_utils.R")

## args = commandArgs(trailingOnly=TRUE)

### SCRIPT FLAGS AND PATH VARIABLES ###

## cal.validation.folder <- args[1]
cal.validation.folder <- "./"
## rdata.path <- args[2]
rdata.path <- "OUMVNLikelihoodTwoOptFBDOneTrait_nonultra.RData"
## n.sim <- as.numeric(args[3])
n.sim <- 100
## job.prefix <- args[4]
job.prefix <- "OUMVNTwoOptFBD"
## n.param <- as.numeric(args[5])
n.param <- 5

## param.names <- strsplit(args[6], ",")[[1]]
param.names <- c("sigmasq", "rv", "theta1", "theta2", "alpha")
## beast.param.names <- strsplit(args[7], ",")[[1]]
beast.param.names <- c("OUSigmaSq", "OURootValue", "OUTheta1", "OUTheta2", "OUAlpha")
## param.labs.preparse <- strsplit(args[8], ",")[[1]]
## prior.means.preparse <- strsplit(args[9], ",")[[1]]
## param.labs <- prior.means <- c()
## for (i in 1:n.param) {
##     param.labs[i] = eval(parse(text=param.labs.preparse[i]))
##     prior.means[i] = eval(parse(text=prior.means.preparse[i]))
## }
## mle.param.names <- strsplit(args[11], ",")[[1]]

param.labs <- c(expression(sigma^2), expression(y[0]), expression(theta[1]), expression(theta[2]), expression(alpha))
mle.param.names <- c("sigmasq.mle", "rv.mle", "theta1.mle", "theta2.mle", "alpha.mle")
prior.means <- c(0.233, 0.0, 1.0, 1.0, 1.1633)

tree.type <- "nonultra"

res.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_results/")
res.files <- mixedsort(paste0(res.path,list.files(res.path)))
load(rdata.path)

############# DOING STUFF #############

n.cols = n.param * 3
log.df <- data.frame(matrix(ncol=n.cols, nrow=n.sim))
names(log.df) <- as.vector(outer(c("lower", "upper", "mean"), param.names, paste, sep="."))
cols <- seq(length.out=n.param, by=3) # 3=lower, upper, mean
cols.2.compare <- c(3, 6)
for (i in 1:n.sim) {
    this.sim.df = read.delim(res.files[i], comment.char='#')

    k = 1
    for (j in cols) {
        beast.param.name = as.character(beast.param.names[k])
        ## print(paste("i=", i, " param=", beast.param.name))
        log.df[i,j:(j+2)] = get.95(this.sim.df[102:1001,beast.param.name])
        k = k+1
    }
}

# putting true values and estimated ones together
full.df <- cbind(true.param.df, log.df)
cols.2.flip <- c(3, 8) # 3 and 4 should flip, 8 and 9 should flip
full.df <- flip.trace.var(full.df, cols.2.flip)
plot.hdi <- vector("list", n.param)
for (i in 1:n.param) {
    upper = paste0("upper.", param.names[i])
    lower = paste0("lower.", param.names[i])
    plot.hdi[[i]] = (full.df[,upper] < full.df[,param.names[i]] | full.df[,lower] > full.df[,param.names[i]])
}

# coverage
table((full.df$sigmasq >= full.df$lower.sigmasq) & (full.df$sigmasq <= full.df$upper.sigmasq))
table((full.df$rv >= full.df$lower.rv) & (full.df$rv <= full.df$upper.rv))
table((full.df$theta1 >= full.df$lower.theta1) & (full.df$theta1 <= full.df$upper.theta1))
table((full.df$theta2 >= full.df$lower.theta2) & (full.df$theta2 <= full.df$upper.theta2))
table((full.df$alpha >= full.df$lower.alpha) & (full.df$alpha <= full.df$upper.alpha))

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
    min.y = min(full.df[,paste0("mean.",param.names[i])], na.rm=TRUE)
    max.y = max(full.df[,paste0("mean.",param.names[i])], na.rm=TRUE)
    all.plots[[i]] = get.plot(param.names[i], paste0("mean.",param.names[i]),
                              min.x, max.x, min.x, max.x, x.lab, prior.means[i],
                              full.df, plot.hdi[[i]])
    ## all.plots[[i]] = get.plot(param.names[i], paste0(param.names[i],".mle"),
    ##                           min.x, max.x, min.x, max.x, x.lab, prior.means[i],
    ##                           full.df)
    names(all.plots)[i] = paste0("plot", i)
}
list2env(all.plots, .GlobalEnv) # sending plots in list into environment so I cna use plot_grid

png(paste0(cal.validation.folder, job.prefix, "_", tree.type, "_graphs.png"), height=24, width=10, unit="cm", res=300)
ggarrange(plotlist=all.plots, ncol=1, nrow=3)
dev.off()
