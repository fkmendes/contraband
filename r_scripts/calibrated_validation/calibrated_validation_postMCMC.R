library(ggplot2)
library(gtools)

source("calibrated_validation_utils.R")

n.sim <- 100
n.param <- 2
sigma.rate <- 5

param.names <- c("sigmasq", "mu")
beast.param.names <- c("BMSigmaSq", "BMMean")

## template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_template.xml"
## res.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_results/"
## res.files <- mixedsort(paste0(res.path,list.files(res.path)))
## rdata.path <- gsub("_template.xml", ".RData", template.path)
## load(rdata.path)

## for non-ultrametric analysis (comment out or in w.r.t. to lines above)
template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_nonultra_template.xml"
res.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_nonultra_results/"
res.files <- mixedsort(paste0(res.path,list.files(res.path)))
rdata.path <- gsub("_template.xml", ".RData", template.path)
load(rdata.path)

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
plot = ggplot() + geom_point(data=full.df, mapping=aes(x=sigmasq, y=mean.sigmasq), shape=20) +
    coord_cartesian(ylim=c(min.x, max.x)) +
    xlab(expression(sigma^2)) + ylab("Posterior mean") +
    geom_abline(slope=1, linetype="dotted") +
    geom_abline(slope=0, intercept=(1/sigma.rate), color="blue") +
    theme(
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.line = element_line(),
        axis.ticks = element_line(color="black"),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)
    ) + scale_x_continuous(labels = function(x) round(as.numeric(x), digits=3), limits=c(min.x,max.y))
plot

min.x <- min(full.df$mu)
max.x <- max(full.df$mu)
min.y <- min(full.df$mean.mu)
max.y <- max(full.df$mean.mu)
plot = ggplot() + geom_point(data=full.df, mapping=aes(x=mu, y=mean.mu), shape=20) +
    coord_cartesian(ylim=c(min.x, max.x)) +
    xlab(expression(mu)) + ylab("Posterior mean") +
    geom_abline(slope=1, linetype="dotted") +
    geom_abline(slope=0, intercept=0, color="blue") +
    theme(
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.line = element_line(),
        axis.ticks = element_line(color="black"),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)
    ) + scale_x_continuous(labels = function(x) round(as.numeric(x), digits=3), limits=c(min.x,max.x))
plot
