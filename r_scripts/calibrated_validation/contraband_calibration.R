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
## sigma.meanlog <- 0.0
## sigma.sdlog <- 1.0
sigma.rate <- 5
x0.mean <- 0.0
x0.sd <- 2.0

template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_template.xml"
xmlfolder.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_xmls/"
## res.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_results/"
shell.scripts.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNOneTrait_shellscripts/"
rdata.path <- gsub("_template.xml", ".RData", template.path)

# for server analyses
xml.file.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/BMMVNOneTrait_xmls/"
xml.file.prefix <- "BMMVNLikelihoodOneTrait_fixedtree_"
res.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/BMMVNOneTrait_results/"
jar.path <- "/nesi/project/nesi00390/fkmendes/contraband/contraband.jar"
time.needed <- "00:15:00"
job.prefix <- "BMMVNOneTrait"

############# DOING STUFF #############

## in ape, branch numbering is done in postorder order: to see the order of branches, you can do tr <- reorder(tr, "po"); plor(tr); edgelabels()

set.seed(123)

## simulating birth-death tree
lambda <- rexp(1, rate=80) # lambda from exponential (mean = 1/80)
mu <- rexp(1, rate=100) # mu from exponential (mean = 1/100)
tr <- sim.bd.taxa.age(n.spp, 1, lambda, mu, age=100, mrca=TRUE)[[1]]; write.tree(tr) # mrca=TRUE means the process starts at the root (i.e., no root edge)

## simulating quant trait data sets
## sigmas <- rlnorm(n.sim, meanlog=sigma.meanlog, sdlog=sigma.sdlog); sigmas[1] # 1.201629
sigmas <- rexp(n.sim, rate=sigma.rate); sigmas[1] # 0.02914
x0s <- rnorm(n.sim, mean=x0.mean, sd=x0.sd); x0s[1] # 0.1836
datasets <- vector("list", n.sim) # storing sims
true.param.df <- data.frame(sigmas, x0s)
names(true.param.df) <- c("mu", "sigmasq")
save(true.param.df, file=rdata.path)
    
## for putting on template
traits.4.template <- vector("list", n.sim) 
taxon.strs.4.template <- paste(paste0("<taxon id=\"", tr$tip.label, "\" spec=\"Taxon\"/>"), collapse="\n                  ")

## actually simulating and populating strs for template
if (simulate) {
    for (i in 1:n.sim) {
        datasets[[i]] = fastBM(tr, sig2=sigmas[i], mean=x0s[i])
        traits.4.template[[i]] = paste(paste0(names(datasets[[i]]), "=", datasets[[i]]), collapse=",")
    }
}

## MLE
## res <- mvBM(tr, dataset, model="BM1")
## res$theta
## res$sigma

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
            ## line = gsub("\\[BMSigmaSqPriorMeanHere\\]", format(sigma.meanlog, nsmall=1), line)
            ## line = gsub("\\[BMSigmaSqPriorStdDevHere\\]", format(sigma.sdlog, nsmall=1), line)
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

if (write.shellscripts) {
    for (sim.idx in 1:n.sim) {
        write.shell.script(shell.scripts.path, sim.idx, time.needed, job.prefix, jar.path, xml.file.path, xml.file.prefix)
    }
}
