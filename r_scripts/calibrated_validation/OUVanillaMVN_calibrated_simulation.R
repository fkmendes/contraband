library(TreeSim)
library(mvMORPH)

source("calibrated_validation_utils.R")

### SCRIPT FLAGS AND PATH VARIABLES ###

simulate <- TRUE
write.xmls <- TRUE
write.shellscripts <- TRUE

n.param <- 4
n.sim <- 100
n.spp <- 50
sigma.rate <- 5
rv.mean <- 0.0 # root value (theta0 in mvMORPH, the first element in the theta vector result)
rv.sd <- 2.0 #
th.mean <- 1.0 # theta1
th.sd <- 2.0 #
th2.mean <- 2.0 # theta2
th2.sd <- 2.0 #
th3.mean <- 3.0 # theta3
th3.sd <- 2.0 #
alpha.mean <- alpha.sd <- 1.0

## template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/OUVanillaMVNLikelihoodOneTrait_fixedtree_template.xml"
## xmlfolder.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/OUVanillaMVNOneTrait_xmls/"
## shell.scripts.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/OUVanillaMVNOneTrait_shellscripts/"
## rdata.path <- gsub("_template.xml", ".RData", template.path)

## # for server analyses
## xml.file.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/OUVanillaMVNOneTrait_xmls/"
## xml.file.prefix <- "OUVanillaMVNLikelihoodOneTrait_fixedtree_"
## res.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/OUVanillaMVNOneTrait_results/"
## jar.path <- "/nesi/project/nesi00390/fkmendes/contraband/contraband.jar"
## time.needed <- "01:00:00"
## job.prefix <- "OUVanilla"

############# DOING STUFF #############

## in ape, branch numbering is done in postorder order: to see the order of branches, you can do tr <- reorder(tr, "po"); plor(tr); edgelabels()

set.seed(123)

## simulating birth-death tree
lambda <- rexp(1, rate=80) # lambda from exponential (mean = 1/80)
mu <- rexp(1, rate=100) # mu from exponential (mean = 1/100)
## tr <- sim.bd.taxa.age(n.spp, 1, lambda, mu, age=100, mrca=TRUE)[[1]]; write.tree(tr) # mrca=TRUE means the process starts at the root (i.e., no root edge)
tr <- read.tree(text="(((((t35:0.1,t32:0.1):0.1,t10:0.1):0.1,t18:0.1):0.1,(((t47:0.1,t9:0.1):0.1,(t43:0.1,t38:0.1):0.1):0.1,(((((t20:0.1,t14:0.1):0.1,t19:0.1):0.1,(t24:0.1,(((t50:0.1,t8:0.1):0.1,t25:0.1):0.1,(t12:0.1,t5:0.1):0.1):0.1):0.1):0.1,t37:0.1):0.1,(t42:0.1,(t13:0.1,t41:0.1):0.1):0.1):0.1):0.1):0.1,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0;")
tr2 <- tr
tr2$tip.label <- paste0(tr$tip.label, "[&amp;Regime=0]")
tr2.newick <- gsub("):", ")[&amp;Regime=0]:", write.tree(tr))
tr2.newick <- gsub("-", ";", tr.newick) # gotta do this because weird dashes appear above

## simulating quant trait data sets
sigmas <- rexp(n.sim, rate=sigma.rate);
rvs <- rnorm(n.sim, mean=rv.mean, sd=rv.sd);
ths <- rnorm(n.sim, mean=th.mean, sd=th.sd);
th2s <- rnorm(n.sim, mean=th2.mean, sd=th2.sd);
alphas <- rlnorm(n.sim, mean=alpha.mean, sd=alpha.sd);

datasets <- vector("list", n.sim) # storing sims
mles <- data.frame(matrix(NA,100,n.param))

## for putting on template
traits.4.template <- vector("list", n.sim) 
taxon.strs.4.template <- paste(paste0("<taxon id=\"", tr$tip.label, "\" spec=\"Taxon\"/>"), collapse="\n              ")

## actually simulating and populating strs for template
if (simulate) {
    for (i in 1:n.sim) {
        ## using mvMORPH
        ### two regimes
        this.tr = paintSubTree(tr, node=52, state=2)
        
        ### fig2 (right, Cecile's paper)
        ## this.tr = paintSubTree(tr, node=60, state=ths[i], anc.state=rvs[i])
        ## this.tr = paintSubTree(this.tr, node=77, state=th2s[i], anc.state=rvs[i])
        ## this.tr = paintSubTree(this.tr, node=87, state=rvs[i])

        ### fig 2 (left, Cecile's paper)
        ## this.tr = paintSubTree(tr, node=52, state=ths[i], anc.state=rvs[i])

        datasets[[i]] = mvSIM(this.tr, param=list(sigma=sigmas[i],
                                                  alpha=alphas[i],
                                                  theta=c(rvs[i], ths[i]),
                                             ## theta=c(rvs[i],ths[i],th2s[i],th3s[i]),
                                             root=TRUE), model="OU1")

        print(paste0("Calling mvOU for sim #",i))
        ## mle.res = mvOU(this.tr, datasets[[i]], model="OUM", param=list(vcv="fixedRoot", root=TRUE))
        mle.res = mvOU(tr, datasets[[i]], model="OU1", param=list(vcv="randomRoot", root=TRUE))
        mles[i,] = c(mle.res$sigma, mle.res$theta, mle.res$alpha)
        traits.4.template[[i]] = paste(paste0(row.names(datasets[[i]]), "=", datasets[[i]]), collapse=",")
    }
}

## saving true values and MLEs to .RData
## true.param.df <- cbind(data.frame(sigmas, rvs, alphas), mles) # theta0 all over
true.param.df <- cbind(data.frame(sigmas, rvs, ths, alphas), mles) # theta0, then theta1 all over;
## true.param.df <- cbind(data.frame(sigmas, rvs, ths, th2s, alphas), mles) # theta0 = old regime, theta1, theta2
## true.param.df <- cbind(data.frame(sigmas, rvs, ths, th2s, th3s, alphas), mles) # theta0, theta1 (old regime), theta2, theta3
## names(true.param.df) <- c("sigmasq", "rv", "alpha", "sigmasq.mle", "rv.mle", "alpha.mle")
names(true.param.df) <- c("sigmasq", "rv", "theta", "alpha", "sigmasq.mle", "rv.mle", "theta.mle", "alpha.mle")
## names(true.param.df) <- c("sigmasq", "rv", "theta1", "theta2", "alpha", "sigmasq.mle", "rv.mle", "theta1.mle", "theta2.mle", "alpha.mle")
## names(true.param.df) <- c("sigmasq", "rv", "theta1", "theta2", "theta3", "alpha", "sigmasq.mle", "rv.mle", "theta1.mle", "theta2.mle", "alpha.mle")
## save(true.param.df, file=rdata.path)

## MLE correlation plot (just for sanity check)
## png("~/Desktop/ultrametric.png", width=15, height=20, unit="cm", res=300)
## png("~/Desktop/non-ultrametric.png", width=15, height=20, unit="cm", res=300)
par(mfrow=c(3,2))
plot(rv.mle~rv, data=true.param.df, xlab="Root value", ylab="MLE of root value", pch=20, ylim=c(min(true.param.df$rv),max(true.param.df$rv)))
## plot(theta.mle~theta, data=true.param.df, xlab=expression(theta), ylab=expression(theta[MLE]), pch=20)
plot(theta1.mle~theta1, data=true.param.df, xlab=expression(theta[1]), ylab=expression(paste("MLE ", theta[1])), pch=20)
plot(theta2.mle~theta2, data=true.param.df, xlab=expression(theta[2]), ylab=expression(paste("MLE" ,theta[2])), pch=20, ylim=c(min(true.param.df$theta2), max(true.param.df$theta2)))
plot(alpha.mle~alpha, data=true.param.df, xlab=expression(alpha), ylab=expression(alpha[MLE]), pch=20, ylim=c(min(true.param.df$alpha),max(true.param.df$alpha)))
plot(sigmasq.mle~sigmasq, data=true.param.df, xlab=expression(sigma^2), ylab=expression(sigma[MLE]^2), pch=20, ylim=c(min(true.param.df$sigmasq),max(true.param.df$sigmasq)))
plot(tr, tip=FALSE)
## plotSimmap(this.tr, ftype="off")
## dev.off()

## writing xmls
## if (write.xmls) {
##     for (sim.idx in 1:n.sim) {
##         xml.file.name = basename(gsub("template.xml", paste0(sim.idx, ".xml"), template.path))
##         if (file.exists(paste0(xmlfolder.path, xml.file.name))) {
##             file.remove(paste0(xmlfolder.path, xml.file.name))
##         }
        
##         file.name = basename(gsub(".xml", "", xml.file.name))

##         template.lines = readLines(template.path)
##         for (line in template.lines) {
##             line = gsub("\\[OUSigmaSqPriorMeanHere\\]", format(1/sigma.rate, nsmall=1), line)
##             line = gsub("\\[OURootValuePriorMeanHere\\]", format(rv.mean, nsmall=1), line)
##             line = gsub("\\[OURootValuePriorStdDevHere\\]", format(rv.sd, nsmall=1), line)
##             line = gsub("\\[OUAlphaPriorMeanHere\\]", format(alpha.mean, nsmall=1), line)
##             line = gsub("\\[OUAlphaPriorStdDevHere\\]", format(alpha.sd, nsmall=1), line)
##             line = gsub("\\[OUThetaPriorMeanHere\\]", format(th.mean, nsmall=1), line)
##             line = gsub("\\[OUThetaPriorStdDevHere\\]", format(th.sd, nsmall=1), line)
##             line = gsub("\\[TreeHere\\]", tr2.newick, line)
##             line = gsub("\\[TaxonSetHere\\]", taxon.strs.4.template, line)
##             line = gsub("\\[TraitValuesHere\\]", traits.4.template[[sim.idx]], line)
##             line = gsub("\\[FileNameHere\\]", paste0(res.path, file.name), line)
##             write(line, file=paste0(xmlfolder.path, xml.file.name), append=TRUE)
##         }
##     }
## }

## writing shell scripts
## if (write.shellscripts) {
##     for (sim.idx in 1:n.sim) {
##         write.shell.script(shell.scripts.path, sim.idx, time.needed, job.prefix, jar.path, xml.file.path, xml.file.prefix)
##     }
## }
