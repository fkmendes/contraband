library(TreeSim)
library(phytools)
library(stringr)
source("calibrated_validation_utils.R")

args = commandArgs(trailingOnly=TRUE)

### SCRIPT FLAGS AND PATH VARIABLES ###

simulate <- args[1]
write.xmls <- args[2]
write.shellscripts <- args[3]
cal.validation.folder <- args[4]
n.sim <- as.numeric(args[5])
n.spp <- as.numeric(args[6])
job.prefix <- args[7] # e.g., "BMMVN" or "BMPrune"
time.needed <- args[8] # e.g., "00:15:00"
template.name <- args[9]
tree.type <- args[10]
xmlfolder.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_xmls/")
xml.file.prefix <- args[11]
shell.scripts.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_shellscripts/")
template.path <- paste0(cal.validation.folder, template.name)
rdata.path <- gsub("_template.xml", ".RData", template.path)

# cluster stuff
cluster.validation.folder <- args[12]
xml.file.path <- paste0(cluster.validation.folder, "BMMVNOneTrait_", tree.type, "_xmls/")
res.path <- paste0(cluster.validation.folder, "BMMVNOneTrait_", tree.type, "_results/")
jar.path <- args[13]

n.param <- 2
sigma.rate <- 5 # exponential prior
x0.mean <- 0.0 # normal prior
x0.sd <- 2.0

## template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_template.xml"
## xmlfolder.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_xmls/"
## shell.scripts.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_shellscripts/"
## rdata.path <- gsub("_template.xml", ".RData", template.path)

# for server analyses
## xml.file.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/BMMVNOneTrait_xmls/"
## xml.file.prefix <- "BMMVNLikelihoodOneTrait_fixedtree_"
## res.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/BMMVNOneTrait_results/"
## jar.path <- "/nesi/project/nesi00390/fkmendes/contraband/contraband.jar"

############# DOING STUFF #############

## in ape, branch numbering is done in postorder order: to see the order of branches, you can do tr <- reorder(tr, "po"); plor(tr); edgelabels()

set.seed(123)

## simulating birth-death tree
lambda <- rexp(1, rate=80) # lambda from exponential (mean = 1/80)
mu <- rexp(1, rate=100) # mu from exponential (mean = 1/100)
## tr <- sim.bd.taxa.age(n.spp, 1, lambda, mu, age=100, mrca=TRUE)[[1]]; # write.tree(tr) # mrca=TRUE means the process starts at the root (i.e., no root edge)
tr <- read.tree(text="(((((t35:0.1,t32:0.1):0.1,t10:0.1):0.1,t18:0.1):0.1,(((t47:0.1,t9:0.1):0.1,(t43:0.1,t38:0.1):0.1):0.1,(((((t20:0.1,t14:0.1):0.1,t19:0.1):0.1,(t24:0.1,(((t50:0.1,t8:0.1):0.1,t25:0.1):0.1,(t12:0.1,t5:0.1):0.1):0.1):0.1):0.1,t37:0.1):0.1,(t42:0.1,(t13:0.1,t41:0.1):0.1):0.1):0.1):0.1):0.1,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0;")

## simulating quant trait data sets
sigmas <- rexp(n.sim, rate=sigma.rate); # sigmas[1] # 0.02914
x0s <- rnorm(n.sim, mean=x0.mean, sd=x0.sd); # x0s[1] # 2.219841
datasets <- vector("list", n.sim) # storing sims
mles <- data.frame(matrix(NA,100,n.param))

## for putting on template
traits.4.template <- vector("list", n.sim)
spnames.4.template <- paste(tr$tip.label, collapse=",")
taxon.strs.4.template <- paste(paste0("<taxon id=\"", tr$tip.label, "\" spec=\"Taxon\"/>"), collapse="\n                  ")

## actually simulating and populating strs for template
if (simulate) {
    for (i in 1:n.sim) {
        datasets[[i]] = fastBM(tr, sig2=sigmas[i], a=x0s[i])
        cat(paste0("Calling fitContinuous for sim #",i,"\n"))
        mle.res = fitContinuous(tr, datasets[[i]], model="BM")
        mles[i,] = c(mle.res$opt$sigsq, mle.res$opt$z0)
        traits.4.template[[i]] = paste(datasets[[i]], collapse=" ")
    }

    ## saving true values and MLEs to .RData
    true.param.df <- cbind(data.frame(sigmas, x0s), mles)
    names(true.param.df) <- c("sigmasq", "mu", "sigmasq.mle", "mu.mle")
    save(tr, true.param.df, file=rdata.path)
} else {
    load(rdata.path) # don't simulate, just load saved simulation
}

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
            line = gsub("\\[SpNamesHere\\]", spnames.4.template, line)
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
