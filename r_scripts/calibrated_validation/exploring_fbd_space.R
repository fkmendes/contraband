library(phytools)
library(FossilSim)

find.success <- function(list.of.trees) {
    trees = 0
    n.sa = c()
    n.extant = c()
    n.tips = c()
    for (i in 1:length(list.of.trees)) {
        if (length(class(list.of.trees[[i]])) != 1) {
            trees = trees+1
            sas = sum(list.of.trees[[i]]$edge.length==0)
            tips = length(list.of.trees[[i]]$tip.label)
            n.sa = c(n.sa, sas)
            n.extant = c(n.extant, tips-sas)
            n.tips = c(n.tips, tips)
        }
    }
    res = vector("list", 4)
    res[[1]] = trees
    res[[2]] = n.sa
    res[[3]] = n.extant
    res[[4]] = n.tips
    return(res)
}

tmp <- sim.fbd.age(1, 100, 3.5, 0.1, .5, complete=TRUE, mrca=TRUE)
find.success(tmp)
