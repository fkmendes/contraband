# author: Fabio K. Mendes
# This R script gives us the expected values for JUnit tests
# (1) 'MVNShiftVcvMatTest'

library(mvMORPH)
library(phytools)

## (1) 'MVNShiftVcvMatTest'

tr <- read.tree(text="((sp1:1.0,sp2:1.0):1.0,sp3:2.0);")
tr <- paintSubTree(tr, node=5, state=2)
plotSimmap(tr)
# black = 0.1, red = 0.2

## sigmas <- matrix(c((0.1+5.0), 0.1, 0.0, 0.1, (0.1+5.0), 0.0, 0.0, 0.0, 0.1*2), 3) # works as well, but not the "right" way
sigmas <- list(black=2.0, red=0.2)

set.seed(1)
dat <- mvSIM(tr, nsim=1, model="BMM", param=list(ntraits=1, sigma=sigmas, theta=0.0))
res <- mvBM(tr, dat, model="BMM")
res$sigma # 0.2804859, 0.226648
res$theta # 0.362281
res$LogLik # -3.106221

## (2) 'MVNShiftVcvMatTest2'

tr <- read.tree(text="(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):2.0,(sp4:2.5,sp5:2.5):1.5);")
## plotTree(tr, branch.numbers=T)
## edgelabels()
## painting internal nodes
tr <- paintSubTree(tr, node=7, state=3, stem=TRUE)
tr <- paintSubTree(tr, node=8, state=2, stem=TRUE)
tr <- paintSubTree(tr, node=9, state=5, stem=TRUE)

## painting tips
tr <- paintBranches(tr, edge=1, state=4)
tr <- paintBranches(tr, edge=2, state=5)
tr <- paintBranches(tr, edge=3, state=1)
tr <- paintBranches(tr, edge=4, state=1)
tr <- paintBranches(tr, edge=5, state=1)
plotSimmap(tr)

## sigmas <- matrix(c(2.2, 1.4, 0.8, 0.0, 0.0, 1.4, 2.2, 0.8, 0.0, 0.0, 0.8, 0.8, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.5, 0.0, 0.0, 0.0, 1.5, 2.0),5)
sigmas <- list()

set.seed(2)
dat <- mvSIM(tr, model="BMM", nsim=1, param=list(ntraits=1, sigma=sigmas, theta=0.0))
res <- mvBM(tr, dat, model="BMM")
res$sigma # 0.05573311, 4.31897, 0.2748663, 1.884965e-13, 0.3798448
res$theta # -1.047182 (root value)
res$LogLik # -6.355337
