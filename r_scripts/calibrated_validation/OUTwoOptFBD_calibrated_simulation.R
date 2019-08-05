library(mvMORPH)
library(FossilSim)
library(phytools)
library(stringr)
library(phangorn)
library(cgwtools) # for resave

source("calibrated_validation_utils.R")

args = commandArgs(trailingOnly=TRUE)
print(args)

### SCRIPT FLAGS AND PATH VARIABLES ###

## simulate <- TRUE
## write.xmls <- TRUE
## write.shellscripts <- TRUE
## cal.validation.folder <- "/home/fkur465/Documents/uoa/contraband/r_scripts/calibrated_validation/"
## n.sim <- 100
## n.spp <- 50
## job.prefix <- "OUMVNTwoOptFBD"
## time.needed <- "24:00:00"
## template.name <- "OUMVNLikelihoodTwoOptFBDOneTrait_nonultra_template.xml"
## tree.type <- "nonultra"

simulate <- args[1]
write.xmls <- args[2]
write.shellscripts <- args[3]
cal.validation.folder <- args[4]
n.sim <- as.numeric(args[5])
n.spp <- as.numeric(args[6])
job.prefix <- args[7]
time.needed <- args[8]
template.name <- args[9]
tree.type <- args[10]

xml.file.prefix <- paste0(job.prefix, "OneTrait_", tree.type, "_")
xmlfolder.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_xmls/")
shell.scripts.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_shellscripts/")
template.path <- paste0(cal.validation.folder, template.name)
rdata.path <- gsub("_template.xml", ".RData", template.path)

# cluster stuff
## cluster.validation.folder <- cal.validation.folder
cluster.validation.folder <- args[11]

xml.file.path <- paste0(cluster.validation.folder, job.prefix, "OneTrait_", tree.type, "_xmls/")
res.path <- paste0(cluster.validation.folder, job.prefix, "OneTrait_", tree.type, "_results/")
## jar.path <- paste0(cluster.validation.folder, "contraband.jar")
jar.path <- args[12]

n.param <- 5

## bm-like
## hpd <- c(0.5, 5.0) # sigsq
## optim(c(0, 1), function(x) sum(qlnorm(p=c(0.025,0.975),meanlog=x[1], sdlog=x[2])-hpd)^2) # 0.0932 0.8005
## hpd <- c(0.001, 0.01) # alpha
## optim(c(0, 1), function(x) sum(qlnorm(p=c(0.025,0.975),meanlog=x[1], sdlog=x[2])-hpd)^2) # -5.9691 0.7171

## ou-like
## sigsq: -5.9691 0.7171 # mean will be 0.003319
## alpha: 0.0932 0.8005 # mean will be 1.511

## sigma.rate <- 5 # exponential prior
sigma.mean <- -5.9691 # ou-like
sigma.sd <- 0.7171
## sigma.mean <- 0.0932 # bm-like
## sigma.sd <- 0.8005
rv.mean <- 0.0 # root value (theta0 in mvMORPH, the first element in the theta vector result)
rv.sd <- 2.0 #
th.mean <- 1.0 # thetas
th.sd <- 2.0 #
## alpha.mean <- -5.9691 # bm-like
## alpha.sd <- 0.7171
alpha.mean <- 0.0932 # ou-like
alpha.sd <- 0.8005
## alpha.mean <- alpha.sd <- 1.0

############# DOING STUFF #############

## in ape, branch numbering is done in postorder order: to see the order of branches, you can do tr <- reorder(tr, "po"); plor(tr); edgelabels()

set.seed(123)

## simulating fossilized birth-death trees
lambda <- rexp(1000, rate=80) # lambda from exponential (mean = 1/80)
mu <- rexp(1000, rate=100) # mu from exponential (mean = 1/100)
psi <- rexp(1000, rate=150) # psi from exponential (mean = 1/150)
trs <- vector("list", n.sim + 50)
trs.heights <- vector("list", n.sim)
trs.edge.mean.lengths <- rep(0, (n.sim + 50))
trs.colors.ns <- rep(0, (n.sim + 50))
trs.ntips <- rep(0, (n.sim + 50))
trs.nedges <- rep(0, (n.sim + 50))
trs.nsas <- rep(0, (n.sim + 50))
fossil.entries <- vector("list", n.sim)

success <- 1
counter <- 1
## add 50 trees just in case some trees have SA at the root (in which case likelihoods cant be computed)
if (simulate) {
    while (success <= n.sim + 50) {
        if (lambda[counter] > mu[counter]) {
            if (counter == 279) { set.seed(234) } # this tree takes up all memory... let's ignore it
            tr = sim.fbd.taxa(n.spp, 1, lambda[counter], mu[counter], psi[counter], complete=TRUE)[[1]]
            n.tips = length(tr$tip.label)

            if (n.tips <= 100) {
                cat(paste(c("Simulating tree", success, "with FBD.\n"), sep=" "))
                sa.idxs = tr$edge[,2][tr$edge.length==0] # get the idx of the daughter end of a branch of length 0 (i.e., a sample ancestor idx)
                root.idx = length(tr$tip.label)+1
                ignore.idx = c(sa.idxs, root.idx)
                ## int.node.idxs = (max(which(node.depth.edgelength(tr)==0))+1):(tr$Nnode + length(tr$tip.label))
                node.idxs.to.sample.from = c(1:(n.tips + tr$Nnode))[-ignore.idx]
                random.int.node = sample(node.idxs.to.sample.from, 1)

                # paint the tree
                # the random local clock model "paints" the branch subtending the node with an indicator shift = true ("1")
                trs[[success]] = paintSubTree(tr, node=random.int.node, state=2, stem=FALSE)

                # collecting tree stats (for all n.sim + 50, later we need to get just the ones used in trait simulation)
                trs.ntips[success] = n.tips
                trs.nedges[success] = nrow(tr$edge)
                trs.nsas[success] = sum(tr$edge.length==0)

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

    trs.colors.ns = as.vector(trs.colors.ns)
    theta.assignments = unlist(as.vector((lapply(trs, function(x) paste(rep(0, 2*length(x$tip.label)-1), collapse=" ")))))
    shift.assignments = unlist(as.vector((lapply(trs, function(x) paste(rep("false", 2*length(x$tip.label)-3), collapse=" ")))))
    shift.assignments = paste0("true ", shift.assignments)
    spnames.4.template = unlist(as.vector((lapply(trs, function(x) paste(x$tip.label, collapse=",")))))
}

## simulating quant trait data sets
set.seed(123)
## sigmas <- rexp(n.sim, rate=sigma.rate);
sigmas <- rlnorm(n.sim, mean=sigma.mean, sd=sigma.sd)
rvs <- rnorm(n.sim, mean=rv.mean, sd=rv.sd);
ths1 <- rnorm(n.sim, mean=th.mean, sd=th.sd);
ths2 <- rnorm(n.sim, mean=th.mean, sd=th.sd);
alphas <- rlnorm(n.sim, mean=alpha.mean, sd=alpha.sd);
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
        datasets[[success]] = mvSIM(trs[[counter]], nsim=1, param=list(sigma=sigmas[success],
                                                                      alpha=alphas[success],
                                                                      theta=c(rvs[success], ths1[success], ths2[success]), root=TRUE), model="OUM")
        cat(paste0("Calling mvOU for sim ",success,", tree=",counter,"\n"))
        mle.res = mvOU(trs[[counter]], datasets[[success]], model="OUM", param=list(vcv="fixedRoot", root=TRUE))
        mles[success,] = c(mle.res$sigma, mle.res$theta[[1]], mle.res$theta[[2]], mle.res$theta[[3]], mle.res$alpha)
        traits.4.template[[success]] = paste(datasets[[success]], collapse=" ")

        ## recording success indexes
        successes[success] = counter
        success = success + 1
        counter = counter + 1
    }

    ## saving true values and MLEs to .RData
    true.param.df = cbind(data.frame(sigmas, rvs, ths1, ths2, alphas), mles)
    names(true.param.df) = c("sigmasq", "rv", "theta1", "theta2", "alpha", "sigmasq.mle", "rv.mle", "theta1.mle", "theta2.mle", "alpha.mle")
    trs.ntips.2.save = trs.ntips[successes]
    trs.heights.2.save = trs.heights[successes]
    trs.edge.mean.lengths.2.save = trs.edge.mean.lengths[successes]
    trs.colors.ns.2.save = trs.colors.ns[successes]
    trs.nedges.2.save = trs.nedges[successes]
    trs.nsas.2.save = trs.nsas[successes]
    trees.2.save = trs[successes]
    save(trees.2.save,
         true.param.df,
         trs.ntips.2.save,
         trs.heights.2.save,
         trs.edge.mean.lengths.2.save,
         trs.colors.ns.2.save,
         trs.nedges.2.save,
         trs.nsas.2.save,
         datasets,
         spnames.4.template,
         traits.4.template,
         theta.assignments,
         shift.assignments,
         file=rdata.path)
} else {
    load(rdata.path) # don't simulate, just load saved simulation
}

## MLE correlation plot (just for sanity check)
## plot(sigmasq.mle~sigmasq, data=true.param.df, xlab=expression(sigma^2), ylab=expression(paste(sigma^2)), pch=20, xlim=c(0, max(true.param.df$sigmasq)), ylim=c(0, max(true.param.df$sigmasq)))
## plot(rv.mle~rv, data=true.param.df, xlab=expression(y[0]), ylab=expression(paste("MLE ", y[0])), pch=20, xlim=c(min(true.param.df$rv), max(true.param.df$rv)), ylim=c(min(true.param.df$rv), max(true.param.df$rv)))
## plot(theta1.mle~theta1, data=true.param.df, xlab=expression(theta[1]), ylab=expression(paste("MLE ", theta[1])), pch=20)
## plot(theta2.mle~theta2, data=true.param.df, xlab=expression(theta[2]), ylab=expression(paste("MLE ", theta[2])), pch=20)
## plot(alpha.mle~alpha, data=true.param.df, xlab=expression(alpha), ylab=expression(paste("MLE ", alpha)), pch=20, xlim=c(0, max(true.param.df$alpha)), ylim=c(0, max(true.param.df$alpha)))

mean.trs.h <- mean(unlist(trs.heights)[successes])
mean.trs.edge.length <- mean(unlist(trs.edge.mean.lengths)[successes]) # we want a 1% genetic distance d (=rt), where t is the mean edge length,
                                        # so the nuc evol rate should be 0.01/mean.trs.edge.length

## simulating seqs
seqs.4.template <- vector("list", n.sim)
nuc.rate <- 0.01/mean.trs.edge.length

if (simulate) {
    for (i in 1:n.sim) {
        seq.datasets[[i]] = phyDat2alignment(simSeq(trees.2.save[[i]], l=2000, type="DNA", rate=nuc.rate))

        seq.string = ""
        sp.count = 1
        for (sp.name in trees.2.save[[i]]$tip.label) {
            seq.string = paste(seq.string, paste0("<sequence id=\"", sp.count, "\" taxon=\"", sp.name, "\" totalcount=\"4\"\nvalue=\"", seq.datasets[[i]]$seq[seq.datasets[[i]]$nam==sp.name], "\"/>\n", sep="\n"))
            sp.count = sp.count + 1
        }
        seqs.4.template[[i]] = seq.string
    }

    resave(seq.datasets, seqs.4.template, nuc.rate, file=rdata.path)
}

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
            ## line = gsub("\\[OUSigmaSqPriorMeanHere\\]", format(1/sigma.rate, nsmall=1), line)
            line = gsub("\\[OUSigmaSqPriorMeanHere\\]", format(sigma.mean, nsmall=1), line)
            line = gsub("\\[OUSigmaSqPriorStdDevHere\\]", format(sigma.sd, nsmall=1), line)
            line = gsub("\\[OURootValuePriorMeanHere\\]", format(rv.mean, nsmall=1), line)
            line = gsub("\\[OURootValuePriorStdDevHere\\]", format(rv.sd, nsmall=1), line)
            line = gsub("\\[OUAlphaPriorMeanHere\\]", format(alpha.mean, nsmall=1), line)
            line = gsub("\\[OUAlphaPriorStdDevHere\\]", format(alpha.sd, nsmall=1), line)
            line = gsub("\\[OUThetaPriorMeanHere\\]", format(th.mean, nsmall=1), line)
            line = gsub("\\[OUThetaPriorStdDevHere\\]", format(th.sd, nsmall=1), line)
            ## line = gsub("\\[ThetaAssignmentsHere\\]", theta.assignments[successes[sim.idx]], line)
            line = gsub("\\[ShiftIndicatorsHere\\]", shift.assignments[successes[sim.idx]], line)
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
