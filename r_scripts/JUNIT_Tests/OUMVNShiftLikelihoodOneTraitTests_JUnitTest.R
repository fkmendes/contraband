# author: Fabio K. Mendes
# This R script gives us the expected values for JUnit tests
# (1) 'OUMVNShiftLikelihoodOneTraitTest'

library(mvMORPH)
library(phytools)

tr <- read.tree(text="(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);")
tr <- paintSubTree(tr, node=7, state=2)
plotSimmap(tr)

sigmas <- list(black=2.0, red=0.2)
thetas <- list(black=0.0, red=0.1)

set.seed(1)
dat <- mvSIM(tr, nsim=1, model="OUM", param=list(ntraits=1, sigma=sigmas, theta=thetas, alpha=0.5, vcv="fixedRoot"))
