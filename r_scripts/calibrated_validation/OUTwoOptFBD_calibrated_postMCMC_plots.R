library(ggplot2)
library(gtools)
library(ggpubr)
library(ape)
library(rowr)
source("calibrated_validation_utils.R")

args = commandArgs(trailingOnly=TRUE)
print(args)

### SCRIPT FLAGS AND PATH VARIABLES ###

tree.type <- "nonultra"

cal.validation.folder <- "./"
rdata.path <- "OUMVNLikelihoodTwoOptFBDOneTrait_nonultra.RData"
n.sim <- 100
job.prefix <- "OUMVNTwoOptFBD"
n.param <- 5
param.names <- c("sigmasq", "rv", "theta1", "theta2", "alpha")
## beast.param.names <- c("OUSigmaSq", "OURootValue", "OUTheta1", "OUTheta2", "OUAlpha")
beast.param.names <- c("OUSigmaSq", "OURootValue", "OUThetaSmall", "OUThetaLarge", "OUAlpha")
param.labs <- c(expression(sigma^2), expression(y[0]), expression(theta[1]), expression(theta[2]), expression(alpha))
mle.param.names <- c("sigmasq.mle", "rv.mle", "theta1.mle", "theta2.mle", "alpha.mle")
## prior.means <- c(0.003297929, 0.0, 1.0, 1.0, 1.504103) ## ou-like (alpha high)
prior.means <- c(1.504103, 0.0, 1.0, 1.0, 0.003297929) ## bm-like (alpha low)
## suffix <- "alphahigh" # ou-like
suffix <- "alphalow" # bm-like

cal.validation.folder <- args[1]
rdata.path <- args[2]
n.sim <- as.numeric(args[3])
job.prefix <- args[4]
n.param <- as.numeric(args[5])
param.names <- strsplit(args[6], ",")[[1]]
beast.param.names <- strsplit(args[7], ",")[[1]]
param.labs.preparse <- strsplit(args[8], ",")[[1]]
prior.means.preparse <- strsplit(args[9], ",")[[1]]
param.labs <- prior.means <- c()
for (i in 1:n.param) {
    param.labs[i] = eval(parse(text=param.labs.preparse[i]))
    prior.means[i] = eval(parse(text=prior.means.preparse[i]))
}
mle.param.names <- strsplit(args[10], ",")[[1]]
suffix <- args[11]

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
    names(this.sim.rate.df) = beast.param.names[3:4]
    this.sim.df = read.delim(res.files[i], comment.char='#')
    this.sim.df = cbind.fill(this.sim.rate.df, this.sim.df)

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
bool.vec <- full.df[,3] > full.df[,4]
cols <- c(17,18,19)
for (i in cols) {
    full.df <- flip.w.bool(full.df, bool.vec, i, i+3)
}

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

## ou like
## sigmasq
## FALSE  TRUE
##    24    76
## rv
## FALSE  TRUE
##    10    90
## theta1
## FALSE  TRUE
##     8    92
## theta2
## FALSE  TRUE
##     9    91
## alpha
## FALSE  TRUE
##    26    74

## bm like
## sigmasq
## FALSE  TRUE
##    2     98
## rv
## FALSE  TRUE
##     5    95
## theta1
## FALSE  TRUE
##    10    90
## theta2
## FALSE  TRUE
##     5    95
## alpha
## FALSE  TRUE
##     5    95

### PLOTS ###
# tree height
hs <- data.frame(unlist(trs.heights.2.save))
names(hs) <- c("heights")

plot.heights <- ggplot(data=hs, aes(x=heights)) + geom_histogram(binwidth=50) + theme_classic() + ylab("Count") + xlab("Tree height") + theme(axis.text.x = element_text(color="black", size=10), axis.text.y = element_text(color="black", size=10), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_heights.pdf"), height=3, width=4)
plot.heights
dev.off()

# extant tips
ts <- data.frame(unlist(trs.ntips.2.save) - trs.nsas.2.save/trs.ntips.2.save)
names(ts) <- c("tips")

plot.leaves <- ggplot(data=ts, aes(x="", y=tips)) + geom_boxplot(notch=TRUE, fill="black") + theme_classic() + ylab("Number of leaves") + xlab("") + theme(axis.ticks.y=element_blank(), axis.text.x = element_text(color="black", size=8)) + coord_flip()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_nleaves.pdf"), height=1.25, width=3)
plot.leaves
dev.off()

# pctg derived regime
pctgs <- data.frame(trs.colors.ns.2.save/trs.nedges.2.save)
names(pctgs) <- c("reg2")
plot.regime <- ggplot(data=pctgs, aes(x="", y=reg2)) + geom_boxplot(notch=TRUE, fill="lightgray") + theme_classic() + ylab("Branches with derived regime (%)") + xlab("") + theme(axis.ticks.y=element_blank(), axis.text.x = element_text(color="black", size=8)) + coord_flip()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_regimes.pdf"), height=1.25, width=3)
plot.regime
dev.off()

# pctg of SA's in all trait value data points
sa.pctgs <- data.frame(trs.nsas.2.save/trs.ntips.2.save)
names(sa.pctgs) <- c("sa")
plot.sa <- ggplot(data=sa.pctgs, aes(x="", y=sa)) + geom_boxplot(notch=TRUE) + theme_classic() + ylab("Sampled ancestor data points (%)") + xlab("") + theme(axis.ticks.y=element_blank(), axis.text.x = element_text(color="black", size=8)) + coord_flip()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_sas.pdf"), height=1.25, width=3)
plot.sa
dev.off()

# trait params
all.plots <- vector("list", n.param)
all.plots.mle <- vector("list", n.param)
for (i in 1:n.param) {
    x.lab = param.labs[i]
    min.x = min(full.df[,param.names[i]])
    max.x = max(full.df[,param.names[i]])
    ## min.y = min(full.df[,paste0("mean.",param.names[i])], na.rm=TRUE)
    ## max.y = max(full.df[,paste0("mean.",param.names[i])], na.rm=TRUE)
    min.y = min(full.df[,paste0("lower.",param.names[i])])
    max.y = max(full.df[,paste0("upper.",param.names[i])])
    all.plots.mle[[i]] = get.plot.no.hdi(param.names[i], paste0(param.names[i],".mle"),
                              min.x, max.x, min.y, max.y, x.lab, prior.means[i],
                              full.df)

    names(all.plots.mle)[i] = paste0("plot.mle", i)

    all.plots[[i]] = get.plot(param.names[i], paste0("mean.",param.names[i]),
                              min.x, max.x, min.y, max.y, x.lab, prior.means[i],
                              full.df, plot.hdi[[i]])

    names(all.plots)[i] = paste0("plot", i)
}
list2env(all.plots.mle, .GlobalEnv) # sending plots in list into environment
list2env(all.plots, .GlobalEnv) # sending plots in list into environment

## png(paste0(cal.validation.folder, job.prefix, "_", tree.type, "_graphs.png"), height=24, width=10, unit="cm", res=300)
## ggarrange(plotlist=all.plots, ncol=1, nrow=3)
## dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_sigsq_mle_", suffix, ".pdf"), height=2, width=2.5)
plot.mle1
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_sigsq_hdis_", suffix, ".pdf"), height=2, width=2.5)
plot1
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_y0_mle_", suffix, ".pdf"), height=2, width=2.5)
plot.mle2
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_y0_hdis_", suffix, ".pdf"), height=2, width=2.5)
plot2
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_th1_mle_", suffix, ".pdf"), height=2, width=2.5)
plot.mle3
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_th1_hdis_", suffix, ".pdf"), height=2, width=2.5)
plot3
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_th2_mle_", suffix, ".pdf"), height=2, width=2.5)
plot.mle4
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_th2_hdis_", suffix, ".pdf"), height=2, width=2.5)
plot4
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_alpha_mle_", suffix, ".pdf"), height=2, width=2.5)
plot.mle5
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_alpha_hdis_", suffix, ".pdf"), height=2, width=2.5)
plot5
dev.off()

## ou-like
## sigsq: -5.9691 0.7171
## alpha: 0.0932 0.8005

prior.df <- data.frame(x=c(rlnorm(100000, mean=0.0932, sd=0.8005), rlnorm(100000, mean=-5.9691, sd=0.7171)), y=c(rep("a",100000), rep("s",100000))) ## ou-like
## prior.df <- data.frame(x=c(rlnorm(100000, mean=0.0932, sd=0.8005), rlnorm(100000, mean=-5.9691, sd=0.7171)), y=c(rep("s",100000), rep("a",100000))) ## bm-like

## ou-like
plot.prior.s <- ggplot(subset(prior.df, y=="s"), aes(x=x)) + geom_density(fill="#63ace5") + xlim(0,.03) + ylab("Density") + xlab(expression(sigma^2)) + theme_classic()
## bm-like
## plot.prior.s <- ggplot(subset(prior.df, y=="s"), aes(x=x)) + geom_density(fill="#63ace5") + xlim(0,20) + ylab("Density") + xlab(expression(sigma^2)) + theme_classic()
    ## theme(
    ##     panel.grid.minor = element_blank(),
    ##     panel.border = element_blank(),
    ##     panel.background = element_blank(),
    ##     plot.background = element_blank(),
    ##     plot.title = element_text(hjust=0.5),
    ##     axis.line = element_line(),
    ##     axis.ticks = element_line(color="black"),
    ##     axis.text.x = element_text(color="black", size=10),
    ##     axis.text.y = element_text(color="black", size=10),
    ##     axis.title.x = element_text(size=12),
    ##     axis.title.y = element_text(size=12)
    ## )

## ou-like
plot.prior.a <- ggplot(subset(prior.df, y=="a"), aes(x=x)) + geom_density(fill="#f6cd61") + xlim(0,20) + ylab("Density") + xlab(expression(alpha)) + theme_classic()
## bm-like
## plot.prior.a <- ggplot(subset(prior.df, y=="a"), aes(x=x)) + geom_density(fill="#f6cd61") + xlim(0,.03) + ylab("Density") + xlab(expression(alpha)) + theme_classic()
    ## theme(
    ##     panel.grid.minor = element_blank(),
    ##     panel.border = element_blank(),
    ##     panel.background = element_blank(),
    ##     plot.background = element_blank(),
    ##     plot.title = element_text(hjust=0.5),
    ##     axis.line = element_line(),
    ##     axis.ticks = element_line(color="black"),
    ##     axis.text.x = element_text(color="black", size=10),
    ##     axis.text.y = element_text(color="black", size=10),
    ##     axis.title.x = element_text(size=12),
    ##     axis.title.y = element_text(size=12)
    ## )

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_alpha_priors_", suffix, ".pdf"), height=2, width=2.5)
plot.prior.a
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_sigmasq_priors_", suffix, ".pdf"), height=2, width=2.5) ## ou-like
## pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_sigmasq_priors_alphalow.pdf"), height=2, width=2.5) ## bm-like
plot.prior.s
dev.off()

pdf(paste0(cal.validation.folder, "figs/", job.prefix, "_", tree.type, "_sigmabyalpha_alphahigh.pdf"), height=2.5, width=6)
sig.by.alpha.ml <- ggplot(data=full.df, aes(x=alpha.mle, sigmasq.mle)) + geom_point() + theme_classic() + xlim(0,6) + ylim(0,.015) + xlab(expression(paste("MLE ", alpha))) + ylab(expression(paste("MLE ", sigma^2)))
sig.by.alpha.p <- ggplot(data=full.df, aes(x=mean.alpha, mean.sigmasq)) + geom_point() + theme_classic() + xlab(expression(paste("Mean posterior ", alpha))) + ylab(expression(paste("Mean posterior ", sigma^2)))
ggarrange(sig.by.alpha.p, sig.by.alpha.ml, ncol=2)
dev.off()
