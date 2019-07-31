## library(TreeSim)
library(mvMORPH)
library(FossilSim)
library(phytools)
library(stringr)
library(phangorn)
library(cgwtools) # for resave

source("calibrated_validation_utils.R")

args = commandArgs(trailingOnly=TRUE)

### SCRIPT FLAGS AND PATH VARIABLES ###

simulate <- TRUE
write.xmls <- TRUE
write.shellscripts <- TRUE
cal.validation.folder <- "/Users/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/"
n.sim <- 100
n.spp <- 50
## job.prefix <- "BMMVNShiftTwoRatesFBD"
job.prefix <- "BMMVNShiftTwoRatesFBDfixed"
time.needed <- "15:00:00"
template.name <- "BMMVNShiftLikelihoodTwoRatesFBDfixedOneTrait_nonultra_template.xml"
tree.type <- "nonultra"
xmlfolder.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_xmls/")

## simulate <- args[1]
## write.xmls <- args[2]
## write.shellscripts <- args[3]
## cal.validation.folder <- args[4]
## n.sim <- as.numeric(args[5])
## n.spp <- as.numeric(args[6])
## job.prefix <- args[7] # e.g., "BMMVNShiftTwoRatesFBD"
## time.needed <- args[8] # e.g., "04:00:00"
## template.name <- args[9]
## tree.type <- args[10]
## xml.file.prefix <- args[11]

xml.file.prefix <- paste0(job.prefix, "OneTrait_", tree.type, "_")
shell.scripts.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_shellscripts/")
template.path <- paste0(cal.validation.folder, template.name)
rdata.path <- gsub("_template.xml", ".RData", template.path)

# cluster stuff
## cluster.validation.folder <- cal.validation.folder
cluster.validation.folder <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/"

## cluster.validation.folder <- args[12]

xml.file.path <- paste0(cluster.validation.folder, job.prefix, "OneTrait_", tree.type, "_xmls/")
res.path <- paste0(cluster.validation.folder, job.prefix, "OneTrait_", tree.type, "_results/")
jar.path <- paste0(cluster.validation.folder, "contraband.jar")
## jar.path <- args[13]

n.param <- 3
## sigma.rate <- 5 # exponential prior
sigma.mean <- 0.0932 # log-normal prior
sigma.sd <- 0.8005
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
trs.heights <- vector("list", n.sim + 50)
trs.edge.mean.lengths <- rep(0, (n.sim + 50))
trs.color.ns <- rep(0, (n.sim + 50))
trs.ntips <- rep(0, (n.sim + 50))
trs.nedges <- rep(0, (n.sim + 50))
fossil.entries <- vector("list", n.sim)

success <- 1
counter <- 1
## add 50 trees just in case some trees have SA at the root (in which case likelihoods cant be computed)
while (success <= n.sim + 50) {
    if (lambda[counter] > mu[counter]) {
        if (counter==279) { set.seed(234) }
        tr = sim.fbd.taxa(n.spp, 1, lambda[counter], mu[counter], psi[counter], complete=TRUE)[[1]]

        if (length(tr$tip.label) <= 100) {
            cat(paste(c("Simulating tree", success, "with FBD.\n"), sep=" "))
            int.node.idxs = (max(which(node.depth.edgelength(tr)==0))+1):(tr$Nnode + length(tr$tip.label))

            # uncomment these lines if you want larger painted clades (careful because this changes model)
            ## n.descendants = 1
            ## while ((n.descendants < (length(tr$tip.label) * 0.2)) | (n.descendants > (length(tr$tip.label) * 0.8))) {
                random.int.node = sample(int.node.idxs, 1)
            ##     n.descendants = length(getDescendants(tr, random.int.node))
            ## }

            trs[[success]] = paintSubTree(tr, node=random.int.node, state=2, stem=FALSE)
            trs.ntips[success] = length(tr$tip.label)
            trs.nedges[success] = nrow(tr$edge)
            depths = node.depth.edgelength(trs[[success]])
            tr.height = as.numeric(as.character(max(depths)))
            trs.heights[[success]] = tr.height
            trs.edge.mean.lengths[success] = mean(tr$edge.length[tr$edge.length!=0.0])
            trs.colors.ns[[success]] = as.integer(table(names(unlist(trs[[success]]$maps)))[2]) ## count is for edges with state=2

            fossil.idxs = depths[1:length(trs[[success]]$tip.label)] < tr.height
            fossil.labels = trs[[success]]$tip.label[fossil.idxs]
            fossil.depths.backw = tr.height - depths[1:length(trs[[success]]$tip.label)][fossil.idxs]
            ## print(fossil.labels)
            ## print(fossil.depths.backw)
            fossil.entries[[success]] = paste(paste(fossil.labels, fossil.depths.backw, sep="="), collapse=",\n")
            ## print(fossil.entries[[success]])
            success = success + 1
        }

        else {
            cat("Tree was too large.\n")
        }
    }
    else {
        cat("lambda < mu.\n")
    }
    counter = counter + 1
}

trs.colors.ns <- as.vector(trs.colors.ns)

rate.assignments <- unlist(as.vector((lapply(trs, function(x) paste(rep(0, 2*length(x$tip.label)-1), collapse=" ")))))
shift.assignments <- unlist(as.vector((lapply(trs, function(x) paste(rep("false", 2*length(x$tip.label)-3), collapse=" ")))))
shift.assignments <- paste0("true ", shift.assignments)
spnames.4.template <- unlist(as.vector((lapply(trs, function(x) paste(x$tip.label, collapse=",")))))

## having a look at tree heights
## hist(unlist(trs.heights), prob=T)
## lines(density(rexp(10000, 1/mean.trs.h)), col="red")

## simulating quant trait data sets and nuc seqs
set.seed(123)
## sigmas.1 <- rexp(n.sim, rate=sigma.rate) # sigmas.tmp1[1] # 0.1686
## sigmas.2 <- rexp(n.sim, rate=sigma.rate) # sigmas.tmp2[1] # 0.2275
sigmas.1 <- rlnorm(n.sim, mean=sigma.mean, sd=sigma.sd)
sigmas.2 <- rlnorm(n.sim, mean=sigma.mean, sd=sigma.sd)
## twice.larger.smaller.idx <- which(((sigmas.tmp1/sigmas.tmp2)>=2 | (sigmas.tmp2/sigmas.tmp1)>=2))[1:100]
## sigmas.tmp1 <- sigmas.tmp1[twice.larger.smaller.idx]
## sigmas.tmp2 <- sigmas.tmp2[twice.larger.smaller.idx]
## sigmas.1 <- c(sigmas.tmp1[sigmas.tmp1 < sigmas.tmp2], sigmas.tmp2[!sigmas.tmp1 < sigmas.tmp2])
## sigmas.2 <- c(sigmas.tmp2[sigmas.tmp1 < sigmas.tmp2], sigmas.tmp1[!sigmas.tmp1 < sigmas.tmp2])

x0s <- rnorm(n.sim, mean=x0.mean, sd=x0.sd); # x0s[1] # -0.4015
datasets <- vector("list", n.sim) # storing trait sims
seq.datasets <- vector("list", n.sim) # storing nuc sims
mles <- data.frame(matrix(NA,100,n.param))

## for putting on template
taxon.strs.4.template <- vector("list", n.sim + 50)
for (i in 1:(n.sim+50)) {
    taxon.strs.4.template[[i]] = paste(paste0("<taxon id=\"", trs[[i]]$tip.label, "\" spec=\"Taxon\"/>"), collapse="\n                  ")
}
traits.4.template <- vector("list", n.sim)

seqs.4.template <- vector("list", n.sim)

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
        sigmas = list(black=sigmas.1[success], red=sigmas.2[success])
        datasets[[success]] = mvSIM(trs[[counter]], nsim=1, model="BMM", param=list(sigma=sigmas, theta=x0s[success])) # simulating trait data
        cat(paste0("Calling mvBM for sim ",success,", tree=",counter,"\n"))
        mle.res = mvBM(trs[[counter]], datasets[[success]], model="BMM")
        mles[success,] = c(mle.res$sigma[[1]], mle.res$sigma[[2]], mle.res$theta)
        traits.4.template[[success]] = paste(datasets[[success]], collapse=" ")

        ## recording success indexes
        successes[success] = counter
        success = success + 1
        counter = counter + 1
    }

    ## saving true values and MLEs to .RData
    true.param.df = cbind(data.frame(sigmas.1, sigmas.2, x0s), mles)
    trs.ntips.2.save = trs.ntips[successes]
    trs.heights.2.save = trs.heights[successes]
    trs.edge.mean.lengths.2.save = trs.edge.mean.lengths[successes]
    trs.colors.ns.2.save = trs.colors.ns[successes]
    trs.nedges.2.save = trs.nedges[successes]
    names(true.param.df) = c("sigma1",  "sigma2", "mu", "sigma1.mle", "sigma2.mle", "mu.mle")
    trees.2.save = trs[successes]
    save(trees.2.save, true.param.df, trs.heights.2.save, trs.edge.mean.lengths.2.save, trs.colors.ns.2.save, trs.nedges.2.save, datasets, file=rdata.path)
} else {
    load(rdata.path) # don't simulate, just load saved simulation
}

## plotting how much of the tree is red
## ggplot(data=tmp.df, aes(x="", y=pctg)) + geom_boxplot(notch=TRUE) + theme_classic() + ylab("% of edges with derived regime") + xlab("") + theme(axis.ticks.y=element_blank()) + coord_flip()

## mle correlation plot (just for sanity check)
## plot(mu.mle~mu, data=true.param.df, xlab=expression(mu), ylab=expression(paste(mu[MLE])), pch=20)
## plot(sigma1.mle~sigma1, data=true.param.df, xlab=expression(sigma[1]^2), ylab=expression(paste(sigma^2[MLE])), pch=20)
## plot(sigma2.mle~sigma2, data=true.param.df, xlab=expression(sigma[2]^2), ylab=expression(paste(sigma^2[MLE])), pch=20)

mean.trs.h <- mean(unlist(trs.heights)[successes])
mean.trs.edge.length <- mean(unlist(trs.edge.mean.lengths)[successes]) # we want a 1% genetic distance d (=rt), where t is the mean edge length,
# so the nuc evol rate should be 0.01/mean.trs.edge.length

## simulating seqs
nuc.rate <- 0.01/mean.trs.edge.length
for (i in 1:n.sim) {
    idx = successes[i]
    seq.datasets[[success]] = phyDat2alignment(simSeq(trs[[idx]], l=2000, type="DNA", rate=nuc.rate))

    seq.string = ""
        sp.count = 1
        for (sp.name in trs[[counter]]$tip.label) {
            seq.string = paste(seq.string, paste0("<sequence id=\"", sp.count, "\" taxon=\"", sp.name, "\" totalcount=\"4\"\nvalue=\"", seq.datasets[[success]]$seq[seq.datasets[[success]]$nam==sp.name], "\"/>\n", sep="\n"))
            sp.count = sp.count + 1
        }
    seqs.4.template[[success]] = seq.string
}

resave(seq.datasets, trsfile=rdata.path)

## writing xmls
if (write.xmls) {
    for (sim.idx in 1:n.sim) {
        ## xml.file.name = basename(gsub("template.xml", paste0(sim.idx, ".xml"), template.path))
        xml.file.name = paste0(xml.file.prefix, sim.idx, ".xml")
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
            ## line = gsub("\\[RateValuesPriorMeanHere\\]", format(1/sigma.rate, nsmall=1), line)
            line = gsub("\\[RateValuesPriorMeanHere\\]", format(sigma.mean, nsmall=1), line)
            line = gsub("\\[RateValuesPriorStdDevHere\\]", format(sigma.sd, nsmall=1), line)
            ## line = gsub("\\[RateAssignmentsHere\\]", rate.assignments[successes[sim.idx]], line)
            line = gsub("\\[ShiftIndicatorsHere\\]", shift.assignments[successes[sim.idx]], line)
            line = gsub("\\[BMMeanPriorMeanHere\\]", format(x0.mean, nsmall=1), line)
            line = gsub("\\[BMMeanPriorStdDevHere\\]", format(x0.sd, nsmall=1), line)
            line = gsub("\\[TreeHere\\]", write.tree(trs[[successes[sim.idx]]]), line)
            line = gsub("\\[SpNamesHere\\]", spnames.4.template[successes[sim.idx]], line)
            line = gsub("\\[FossilSetHere\\]", fossil.entries[successes[sim.idx]], line)
            line = gsub("\\[TaxonSetHere\\]", taxon.strs.4.template[successes[sim.idx]], line)
            line = gsub("\\[TraitValuesHere\\]", traits.4.template[[sim.idx]], line)
            line = gsub("\\[FileNameHere\\]", paste0(res.path, file.name), line)
            line = gsub("\\[SeqsHere\\]", seqs.4.template[sim.idx], line)
            line = gsub("\\[MutRateHere\\]", nuc.rate, line)
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
