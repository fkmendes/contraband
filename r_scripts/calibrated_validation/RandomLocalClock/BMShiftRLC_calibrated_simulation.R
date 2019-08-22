library(mvMORPH)
library(FossilSim)
library(phytools)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

template.path <- args[1]
xmlfolder.path <- args[2]
load(args[3])
output.data.folder <- args[4]
n.sim <- as.numeric(args[5])
lambda <- as.numeric(args[6])

#template.path <- "~/Desktop/random_local/BM_Calibrated_Template.xml"
#xmlfolder.path <- "~/Desktop/random_local/XMLs/"
#unmappedtree.path <- "~/Desktop/random_local/unmappedtrs.RData"
#output.data.folder <- "~/Desktop/random_local/output/"
#n.sim <- 200
#nshift.lambda <- 5

sigma.mean <- 0.0932
sigma.sd <- 0.8005
x0s <- 0.0 # fixing root value

# in case we want to run it on a few different trees
tree.nr <- c()
for (i in 1:length(trs)){
  tree.nr[i] <- paste0("tr",i)
}

set.seed(123)
for (r.sim in 1:length(trs)) {
   tr = trs[[r.sim]] # we're only looking at a single tree, but this will do multiple trees if necessary (but .RData passed to this script must have multiple trees inside trs)
   tr.nr = tree.nr[r.sim]

   mapped.trs = vector("list", n.sim)
   datasets = vector("list", n.sim)
   nodes.painted = vector("list", n.sim)

   nshift = rpois(n.sim, lambda) # drawing n shifts from Poisson

   # drawing sigmas (!= number in each simulation) from log-normal
   sigmas = vector("list", n.sim)
   for(i in 1:n.sim) {
       sigmas[[i]] = rlnorm(nshift[i]+1, meanlog=sigma.mean, sdlog=sigma.sd)
   }

   n.tips = length(tr$tip.label)
   sa.idxs = tr$edge[,2][tr$edge.length==0] # get the idx of the daughter end of a branch of length 0 (i.e., a sample ancestor idx)
   root.idx = length(tr$tip.label)+1
   ignore.idx = c(sa.idxs, root.idx)
   node.idxs.to.sample.from = c(1:(n.tips + tr$Nnode))[-ignore.idx]

   set.seed(123)

   # simulate data
   success = 1
   for (i in 1:n.sim) {
      k = nshift[i] # number of shifts in this simulation

      # if there's at least 1 shift, we need to paint the tree
      if (k > 0) {
          t.tr = tr
          random.int.node = sample(node.idxs.to.sample.from, k) # randomly picking an internal node
          nodes.painted[[success]] = random.int.node # store the painted node for later

          node.height = matrix() # we need to record node heights so that painting
          # is done correctly (from most inclusive node to least inclusive)
          for (j in 1:k) {
              node.height[j] = nodeheight(tr, random.int.node[j])
          }

          order.int.node = order(node.height) # storing node indexes in order for painting the oldest nodes first
          color = 2 # start color number at 2

          for (idx in order.int.node) {
              t.tr = paintSubTree(t.tr, node=random.int.node[idx], state=color, stem=TRUE)
              color = color + 1
          }
          ## for(node.nr in 1:k) {
          ##     t.tr = paintSubTree(t.tr, node=random.int.node[order.int.node[node.nr]], state=color, stem=TRUE)
          ##     color <- color + 1
          ## }

          datasets[[i]] = mvSIM(t.tr, nsim=1, model="BMM", param=list(sigma=sigmas[[i]], theta=x0s)) # simulate continuous trait
          mapped.trs[[success]] = t.tr # save painted tree
          success = success + 1
      }

      # if no shifts, different model
      else {
          datasets[[i]] = mvSIM(tr, nsim=1, model="BM1", param=list(sigma=sigmas[[i]], theta=x0s))
          mapped.trs[[success]] = tr
          success = success + 1
      }
   }

   # write .xml files
   taxa = paste(paste0("<taxon id=\"", tr$tip.label, "\" spec=\"Taxon\"/>"), collapse="\n                  ")
   spnames = paste(tr$tip.label, collapse=",")
   template.lines = readLines(template.path)

   for (sim.idx in 1:n.sim) {
       out.file.name = paste0("BMRLC_", tr.nr, "_", sim.idx) # file prefix for naming .log and .trees files later
       xml.file.name = paste0("BMRLC_", tr.nr, "_", sim.idx, ".xml")

       if (file.exists(paste0(xmlfolder.path, xml.file.name))) { print("Deleted already existing .xml file."); unlink(paste0(xmlfolder.path, xml.file.name)) }

       for (line in template.lines) {
           line = gsub("\\[TreeHere\\]", write.tree(tr), line)
           line = gsub("\\[SpNamesHere\\]", spnames, line)
           line = gsub("\\[TaxonSetHere\\]", taxa, line)
           line = gsub("\\[TraitValuesHere\\]", paste(datasets[[sim.idx]], collapse=" "), line)
           line = gsub("\\[LambdaHere\\]", lambda, line)
           line = gsub("\\[LogMeanHere\\]", sigma.mean, line)
           line = gsub("\\[LogSdHere\\]", sigma.sd, line)
           line = gsub("\\[Filename\\]", out.file.name, line)
           write(line, file=paste0(xmlfolder.path, xml.file.name), append=TRUE)
       }
   }

   save.image(paste0(output.data.folder, "BMRLC_simdata_", tr.nr, ".RData"))
}
