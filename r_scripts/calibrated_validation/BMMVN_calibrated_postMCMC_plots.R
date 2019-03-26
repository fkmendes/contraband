library(ggplot2)
library(gtools)
library(gridExtra)
library(sjPlot)

source("calibrated_validation_utils.R")
load(rdata.path)

n.sim <- 100
n.param <- 2
sigma.rate <- 5
x0s.mean <- 0

param.labs <- c(expression(sigma^2), expression(mu))
param.names <- c("sigmasq", "mu")
beast.param.names <- c("BMSigmaSq", "BMMean")
mle.param.names <- c("sigmasq.mle", "mu.mle")
prior.means <- c(1/sigma.rate, x0s.mean)

template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_template.xml"
res.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_results/"
res.files <- mixedsort(paste0(res.path,list.files(res.path)))
rdata.path <- gsub("_template.xml", ".RData", template.path)

## for non-ultrametric analysis (comment out or in w.r.t. to lines above)
## template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_nonultra_template.xml"
## res.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_nonultra_results/"
## res.files <- mixedsort(paste0(res.path,list.files(res.path)))

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

## ultrametric
## mu
##   FALSE  TRUE 
##       6    94
## sigma
##   FALSE  TRUE 
##      10    90 

## ultrametric
## mu
##   FALSE  TRUE 
##       5    95
## sigma
##   FALSE  TRUE 
##      10    90 

### PLOTS ###

all.plots <- vector("list", n.param)
for (i in 1:n.param) {
    x.lab = param.labs[i]
    min.x = min(full.df[,param.names[i]])
    max.x = max(full.df[,param.names[i]])
    min.y = min(full.df[,paste0("mean.",param.names[i])])
    max.y = max(full.df[,paste0("mean.",param.names[i])])
    print(param.names[i]);     print(paste0("mean.",param.names[i]))
    all.plots[[i]] = get.plot(param.names[i], paste0("mean.",param.names[i]),
                              min.x, max.x, min.y, max.y, x.lab, prior.means[i],
                              full.df)
    names(all.plots)[i] = paste0("plot", i)
}
list2env(all.plots, .GlobalEnv) # sending plots in list into environment so I cna use plot_grid

png("~/Desktop/ultrametric_tree.png", height=10, width=10, unit="cm", res=100)
plot(tr, show.tip.label=F)
dev.off()

png("~/Desktop/fossil_tree.png", height=10, width=10, unit="cm", res=100)
plot(tr, show.tip.label=F)
dev.off()

png("~/Desktop/ultrametric_tree_BM_validation.png", height=10, width=20, unit="cm", res=300)
plot_grid(plotlist=mget(paste0("plot", 1:2)))
dev.off()

## png("~/Desktop/fossil_tree_BM_validation.png", height=10, width=20, unit="cm", res=300)
## plot_grid(plotlist=mget(paste0("plot", 1:2)))
## dev.off()

