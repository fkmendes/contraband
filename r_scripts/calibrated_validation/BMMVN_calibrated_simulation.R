library(mvMORPH)
library(TreeSim)
library(phytools)
library(stringr)

source("calibrated_validation_utils.R")

### SCRIPT FLAGS AND PATH VARIABLES ###

simulate <- TRUE
write.xmls <- TRUE
write.shellscripts <- TRUE

n.sim <- 100
n.spp <- 50
sigma.rate <- 5
x0.mean <- 0.0
x0.sd <- 2.0

template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_template.xml"
xmlfolder.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_xmls/"
shell.scripts.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_shellscripts/"
rdata.path <- gsub("_template.xml", ".RData", template.path)

## for non-ultrametric analysis (comment out or in w.r.t. to lines above)
## template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_nonultra_template.xml"
## xmlfolder.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_nonultra_xmls/"
## shell.scripts.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_nonultra_shellscripts/"
## rdata.path <- gsub("_template.xml", ".RData", template.path)

# for server analyses
xml.file.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/BMMVNOneTrait_xmls/"
xml.file.prefix <- "BMMVNLikelihoodOneTrait_fixedtree_"
res.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/BMMVNOneTrait_results/"
jar.path <- "/nesi/project/nesi00390/fkmendes/contraband/contraband.jar"
time.needed <- "00:15:00"
job.prefix <- "BMMVNOneTrait"

## for non-ultrametric analysis (comment out or in w.r.t. to lines above)
## xml.file.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/BMMVNOneTrait_nonultra_xmls/"
## xml.file.prefix <- "BMMVNLikelihoodOneTrait_fixedtree_nonultra_"
## res.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/BMMVNOneTrait_nonultra_results/"
## jar.path <- "/nesi/project/nesi00390/fkmendes/contraband/contraband.jar"
## time.needed <- "00:15:00"


############# DOING STUFF #############

## in ape, branch numbering is done in postorder order: to see the order of branches, you can do tr <- reorder(tr, "po"); plor(tr); edgelabels()

set.seed(123)

## simulating birth-death tree
lambda <- rexp(1, rate=80) # lambda from exponential (mean = 1/80)
mu <- rexp(1, rate=100) # mu from exponential (mean = 1/100)
tr <- sim.bd.taxa.age(n.spp, 1, lambda, mu, age=100, mrca=TRUE)[[1]]; write.tree(tr) # mrca=TRUE means the process starts at the root (i.e., no root edge)
## tr <- read.tree(text="(((((t36:0.1,t15:0.1):0.1,((t32:0.1,t38:0.1):0.1,t30:0.1):0.1):0.1,t48:0.1):0.1,(((((t29:0.1,((t13:0.1,t10:0.1):0.1,((t49:0.1,t44:0.1):0.1,(t12:0.1,(t42:0.1,t28:0.1):0.1):0.1):0.1):0.1):0.1,((t14:0.1,t35:01):0.1,t45:0.1):0.1):0.1,(t8:0.1,(t6:0.1,t47:0.1):0.1):0.1):0.1,((t19:0.1,((t50:0.1,t1:0.1):0.1,t41:0.1):0.1):0.1,t25:0.1):0.1):0.1,t7:0.1):0.1):0.1,((((t40:5.547816064,t39:5.547816064):69.34457458,(t5:63.2303403,(t3:8.905705974,(t2:8.507831451,t24:8.507831451):0.3978745225):54.32463433):11.66205034):2.166829643,t43:77.05922028):20.83631922,(((t21:35.66443979,(t31:19.79691334,t22:19.79691334):15.86752645):43.35619033,(t16:0.798839978,t27:0.798839978):78.22179014):7.146714057,(((t17:11.8481771,(t4:8.654399431,(t33:4.946918061,t20:4.946918061):3.70748137):3.193777673):60.81563855,(t37:26.39721604,(t26:19.7151034,t9:19.7151034):6.682112638):46.26659961):2.831035132,((t23:4.660327926,t18:4.660327926):67.53247256,((t46:6.150543246,t11:6.150543246):29.58028484,t34:35.73082809):36.4619724):3.302050294):10.6724934):11.72819533):2.104460492):0;")

## simulating quant trait data sets
## sigmas <- rlnorm(n.sim, meanlog=sigma.meanlog, sdlog=sigma.sdlog); sigmas[1] # 1.201629
sigmas <- rexp(n.sim, rate=sigma.rate); sigmas[1] # 0.02914
x0s <- rnorm(n.sim, mean=x0.mean, sd=x0.sd); x0s[1] # 0.1836
datasets <- vector("list", n.sim) # storing sims
mles <- data.frame(matrix(NA,100,2))

## for putting on template
traits.4.template <- vector("list", n.sim) 
taxon.strs.4.template <- paste(paste0("<taxon id=\"", tr$tip.label, "\" spec=\"Taxon\"/>"), collapse="\n                  ")

## actually simulating and populating strs for template
if (simulate) {
    for (i in 1:n.sim) {
        datasets[[i]] = fastBM(tr, sig2=sigmas[i], a=x0s[i])
        print(paste0("Calling fitContinuous for sim #",i))
        mle.res = fitContinuous(tr, datasets[[i]], model="BM")
        mles[i,] = c(mle.res$opt$sigsq, mle.res$opt$z0)
        traits.4.template[[i]] = paste(paste0(names(datasets[[i]]), "=", datasets[[i]]), collapse=",")
    }
}

## saving true values and MLEs to .RData
true.param.df <- cbind(data.frame(sigmas, x0s), mles)
names(true.param.df) <- c("sigmasq", "mu", "sigmasq.mle", "mu.mle")
save(true.param.df, file=rdata.path)

## MLE correlation plot (just for sanity check)
## par(mfrow=c(1,2))
## plot(mu.mle~mu, data=true.param.df, xlab=expression(mu), ylab=expression(paste(mu[MLE])), pch=20)
## plot(tr)

## par(mfrow=c(1,2))
## plot(sigmasq.mle~sigmasq, data=true.param.df, xlab=expression(sigma^2), ylab=expression(paste(sigma^2[MLE])), pch=20)
## plot(tr)

## writing xmls
if (write.xmls) {
    for (sim.idx in 1:n.sim) {
        xml.file.name = basename(gsub("template.xml", paste0(sim.idx, ".xml"), template.path))
        if (file.exists(paste0(xmlfolder.path, xml.file.name))) {
            file.remove(paste0(xmlfolder.path, xml.file.name))
        }
        
        file.name = basename(gsub(".xml", "", xml.file.name))

        template.lines = readLines(template.path)
        for (line in template.lines) {
            line = gsub("\\[BMSigmaSqPriorMeanHere\\]", format(1/sigma.rate, nsmall=1), line)
            line = gsub("\\[BMMeanPriorMeanHere\\]", format(x0.mean, nsmall=1), line)
            line = gsub("\\[BMMeanPriorStdDevHere\\]", format(x0.sd, nsmall=1), line)
            line = gsub("\\[TreeHere\\]", write.tree(tr), line)
            line = gsub("\\[TaxonSetHere\\]", taxon.strs.4.template, line)
            line = gsub("\\[TraitValuesHere\\]", traits.4.template[[sim.idx]], line)
            line = gsub("\\[FileNameHere\\]", paste0(res.path, file.name), line)
            write(line, file=paste0(xmlfolder.path, xml.file.name), append=TRUE)
        }
    }
}

## writing shell scripts
if (write.shellscripts) {
    for (sim.idx in 1:n.sim) {
        write.shell.script(shell.scripts.path, sim.idx, time.needed, job.prefix, jar.path, xml.file.path, xml.file.prefix)
    }
}
