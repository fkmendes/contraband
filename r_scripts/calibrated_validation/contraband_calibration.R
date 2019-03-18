library(mvMORPH)
library(TreeSim)
library(ape)
library(stringr)

n.sim <- 100
n.spp <- 50
sigma.meanlog <- 0.0
sigma.sdlog <- 1.0
x0.mean <- 0.0
x0.sd <- 1.0

## in ape, branch numbering is done in postorder order: to see the order of branches, you can do tr <- reorder(tr, "po"); plor(tr); edgelabels()

set.seed(123)

## simulating birth-death tree
lambda <- rexp(1, rate=80) # lambda from exponential (mean = 1/80)
mu <- rexp(1, rate=100) # mu from exponential (mean = 1/100)
tr <- sim.bd.taxa.age(n.spp, 1, lambda, mu, age=100, mrca=TRUE)[[1]]; write.tree(tr) # mrca=TRUE means the process starts at the root (i.e., no root edge)

tr$edge.labels
paste0("<taxon id="

## simulating quant trait data sets
sigmas <- rlnorm(n.sim, meanlog=sigma.meanlog, sdlog=sigma.sdlog); sigmas[1] # 1.201629
x0s <- rnorm(n.sim, mean=x0.mean, sd=x0.sd); x0s[1] # 0.1836
datasets <- vector("list", n.sim) # storing sims

## for putting on template
traits.4.template <- vector("list", n.sim) 
taxon.strs.4.template <- paste(paste0("<taxon id=\"", tr$tip.label, "\" spec=\"Taxon\"/>"), collapse="\n                  ")

## actually simulating and populating strs for template
for (i in 1:n.sim) {
    datasets[[i]] = rTraitCont(tr, model="BM", sigma=sigmas[i], theta=x0s[i])
    traits.4.template[[i]] = paste(paste0(names(datasets[[i]]), "=", datasets[[i]]), collapse=",")
}

## MLE
## res <- mvBM(tr, dataset, model="BM1")
## res$theta
## res$sigma

template.path <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_template.xml"
for (sim.idx in 1:2) {
    xml.file.name = gsub("template.xml", paste0(sim.idx, ".xml"), template.path)
    if (file.exists(xml.file.name)) {
        file.remove(xml.file.name)
    }
    
    file.name = gsub(".xml", "", xml.file.name)

    template.lines = readLines(template.path)
    for (line in template.lines) {
        line = gsub("\\[BMSigmaSqPriorMeanHere\\]", format(sigma.meanlog, nsmall=1), line)
        line = gsub("\\[BMSigmaSqPriorStdDevHere\\]", format(sigma.sdlog, nsmall=1), line)
        line = gsub("\\[BMMeanPriorMeanHere\\]", format(x0.mean, nsmall=1), line)
        line = gsub("\\[BMMeanPriorStdDevHere\\]", format(x0.sd, nsmall=1), line)
        line = gsub("\\[TreeHere\\]", write.tree(tr), line)
        line = gsub("\\[TaxonSetHere\\]", taxon.strs.4.template, line)
        line = gsub("\\[TraitValuesHere\\]", traits.4.template[[sim.idx]], line)
        line = gsub("\\[FileNameHere\\]", file.name, line)
        write(line, file=xml.file.name, append=TRUE)
    }
}


template <- readLines("/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_template.xml")
parse.ith.sim.template(1, "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/BMMVNLikelihoodOneTrait_fixedtree_template.xml", "./")
