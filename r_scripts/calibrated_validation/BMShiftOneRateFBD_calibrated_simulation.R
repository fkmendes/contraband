## library(TreeSim)
library(mvMORPH)
library(FossilSim)
library(phytools)
library(stringr)
source("calibrated_validation_utils.R")

args = commandArgs(trailingOnly=TRUE)

### SCRIPT FLAGS AND PATH VARIABLES ###

## simulate <- TRUE
## write.xmls <- TRUE
## write.shellscripts <- TRUE
## cal.validation.folder <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/"
## n.sim <- 100
## n.spp <- 50
## job.prefix <- "BMMVNShiftOneRateFBD"
## time.needed <- "03:00:00"
## template.name <- "BMMVNShiftLikelihoodOneRateFBDOneTrait_nonultra_template.xml"
## tree.type <- "nonultra"

simulate <- args[1]
write.xmls <- args[2]
write.shellscripts <- args[3]
cal.validation.folder <- args[4]
n.sim <- as.numeric(args[5])
n.spp <- as.numeric(args[6])
job.prefix <- args[7] # e.g., "BMMVNShiftOneRate" or "BMPruneShiftOneRate"
time.needed <- args[8] # e.g., "00:15:00"
template.name <- args[9]
tree.type <- args[10]
xmlfolder.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_xmls/")
xml.file.prefix <- args[11]
xml.file.prefix <- paste0("BMMVNShiftLikelihoodOneRateFBDOneTrait_", tree.type, "_")
shell.scripts.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_shellscripts/")
template.path <- paste0(cal.validation.folder, template.name)
rdata.path <- gsub("_template.xml", ".RData", template.path)

# cluster stuff
## cluster.validation.folder <- cal.validation.folder
## cluster.validation.folder <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/"
cluster.validation.folder <- args[12]
xml.file.path <- paste0(cluster.validation.folder, "BMMVNShiftOneRateFBDOneTrait_", tree.type, "_xmls/")
res.path <- paste0(cluster.validation.folder, "BMMVNShiftOneRateFBDOneTrait_", tree.type, "_results/")
jar.path <- paste0(cluster.validation.folder, "contraband.jar")
## jar.path <- args[13]

n.param <- 2
sigma.rate <- 5 # exponential prior
x0.mean <- 0.0 # normal prior
x0.sd <- 2.0

############# DOING STUFF #############

## in ape, branch numbering is done in postorder order: to see the order of branches, you can do tr <- reorder(tr, "po"); plor(tr); edgelabels()

set.seed(123)

## simulating fossilized birth-death trees
lambda <- rexp(1000, rate=80) # lambda from exponential (mean = 1/80)
mu <- rexp(1000, rate=100) # mu from exponential (mean = 1/100)
psi <- rexp(1000, rate=150) # psi from exponential (mean = 1/150)
trs <- vector("list", n.sim + 50)
trs.heights <- vector("list", n.sim)
fossil.entries <- vector("list", n.sim)

success <- 1
counter <- 1
## add 50 trees just in case some trees have SA at the root (in which case likelihoods cant be computed)
while (success <= n.sim + 50) {
    if (lambda[counter] > mu[counter]) {
        if (counter==276) { set.seed(234) }
        tr = sim.fbd.taxa(n.spp, 1, lambda[counter], mu[counter], psi[counter], complete=TRUE)[[1]]

        if (length(tr$tip.label) <= 100) {
            cat(paste(c("Simulating tree", success, "with FBD.\n"), sep=" "))
            trs[[success]] = tr
            ## trs[[success]] = sim.fbd.taxa(n.spp, 1, lambda[counter], mu[counter], psi[counter], complete=TRUE)[[1]]
            depths = node.depth.edgelength(trs[[success]])
            tr.height = max(depths)
            trs.heights[[success]] = tr.height
            fossil.idxs = depths[1:length(trs[[success]]$tip.label)] < as.character(tr.height)
            fossil.labels = trs[[success]]$tip.label[fossil.idxs]
            fossil.depths.backw = tr.height - depths[1:length(trs[[success]]$tip.label)][fossil.idxs]
            ## print(fossil.labels)
            ## print(fossil.depths.backw)
            fossil.entries[[success]] = paste(paste(fossil.labels, fossil.depths.backw, sep="="), collapse=",\n")
            ## print(fossil.entries[[success]])
            success = success + 1
        }

        else {
            print("Tree was too large")
        }
    }
    else {
        print("lambda < mu.")
    }
    counter = counter + 1
}
rate.assignments <- unlist(as.vector((lapply(trs, function(x) paste(rep(0, 2*length(x$tip.label)-1), collapse=" ")))))
spnames.4.template <- unlist(as.vector((lapply(trs, function(x) paste(x$tip.label, collapse=",")))))
mean.trs.h <- mean(unlist(trs.heights))

## having a look at tree heights
## hist(unlist(trs.heights), prob=T)
## lines(density(rexp(10000, 1/mean.trs.h)), col="red")

## simulating quant trait data sets
sigmas <- rexp(n.sim, rate=sigma.rate); # sigmas[1] # 0.02914
x0s <- rnorm(n.sim, mean=x0.mean, sd=x0.sd); # x0s[1] # 2.219841
datasets <- vector("list", n.sim) # storing sims
mles <- data.frame(matrix(NA,100,n.param))

## for putting on template
taxon.strs.4.template <- vector("list", n.sim + 50)
for (i in 1:(n.sim+50)) {
    taxon.strs.4.template[[i]] = paste(paste0("<taxon id=\"", trs[[i]]$tip.label, "\" spec=\"Taxon\"/>"), collapse="\n                  ")
}
traits.4.template <- vector("list", n.sim) 

## actually simulating and populating strs for template
set.seed(123)
counter <- 1
success <- 1
successes <- rep(1, 100) # indexes of trees in trs that could have likelihoods computed
if (simulate) {
    while (success <= 100) {
        ## BM barfs when the first SA is on the root
        if (sum(apply(nodeHeights(trs[[counter]]),1,sum)==0) >=1) {
            counter = counter + 1
            next
        }
        datasets[[success]] = fastBM(trs[[counter]], sig2=sigmas[success], a=x0s[success])
        cat(paste0("Calling mvBM for sim ",success,", tree=",counter,"\n"))
        mle.res = mvBM(trs[[counter]], datasets[[success]], model="BM1")
        mles[success,] = c(mle.res$sigma, mle.res$theta)
        traits.4.template[[success]] = paste(datasets[[success]], collapse=" ")
        successes[success] = counter
        success = success + 1
        counter = counter + 1
    }

    ## saving true values and MLEs to .RData
    true.param.df = cbind(data.frame(sigmas, x0s), mles)
    names(true.param.df) = c("sigmasq", "mu", "sigmasq.mle", "mu.mle")
    trees.2.save = trs[successes]
    save(trees.2.save, true.param.df, file=rdata.path)
} else {
    load(rdata.path) # don't simulate, just load saved simulation
}

## MLE correlation plot (just for sanity check)
## plot(mu.mle~mu, data=true.param.df, xlab=expression(mu), ylab=expression(paste(mu[MLE])), pch=20)
## plot(sigmasq.mle~sigmasq, data=true.param.df, xlab=expression(sigma^2), ylab=expression(paste(sigma^2[MLE])), pch=20)

## writing xmls
if (write.xmls) {
    for (sim.idx in 1:n.sim) {
        xml.file.name = basename(gsub("template.xml", paste0(sim.idx, ".xml"), template.path))
        print (xml.file.name)
        if (file.exists(paste0(xmlfolder.path, xml.file.name))) {
            file.remove(paste0(xmlfolder.path, xml.file.name))
        }

        file.name = basename(gsub(".xml", "", xml.file.name))

        template.lines = readLines(template.path)

        for (line in template.lines) {
            line = gsub("\\[InitialOriginValue\\]", format((trs.heights[[successes[sim.idx]]]+0.01), nsmall=1), line)
            line = gsub("\\[FBDbirthRatePriorMeanHere\\]", format(1/80, nsmall=1), line)
            line = gsub("\\[FBDdeathRatePriorMeanHere\\]", format(1/100, nsmall=1), line)
            line = gsub("\\[FBDsamplingRatePriorMeanHere\\]", format(1/150, nsmall=1), line)
            line = gsub("\\[FBDoriginPriorMeanHere\\]", mean.trs.h, line)
            line = gsub("\\[RateValuesPriorMeanHere\\]", format(1/sigma.rate, nsmall=1), line)
            line = gsub("\\[RateAssignmentsHere\\]", rate.assignments[successes[sim.idx]], line)
            line = gsub("\\[BMMeanPriorMeanHere\\]", format(x0.mean, nsmall=1), line)
            line = gsub("\\[BMMeanPriorStdDevHere\\]", format(x0.sd, nsmall=1), line)
            line = gsub("\\[TreeHere\\]", write.tree(trs[[successes[sim.idx]]]), line)
            line = gsub("\\[SpNamesHere\\]", spnames.4.template[successes[sim.idx]], line)
            line = gsub("\\[FossilSetHere\\]", fossil.entries[successes[sim.idx]], line)
            line = gsub("\\[TaxonSetHere\\]", taxon.strs.4.template[successes[sim.idx]], line)
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
