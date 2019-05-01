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

## (3) 'BMMVNShiftLikelihoodOneTraitTest'
### See expectations from BMMVNLikelihoodOnetraitTest

## (4) 'BMMVNShiftLikelihoodOneTraitTest2'
set.seed(123)
lambda <- rexp(1, rate=80) # lambda from exponential (mean = 1/80)
mu <- rexp(1, rate=100) # mu from exponential (mean = 1/100)
tr <- sim.bd.taxa.age(n.spp, 1, lambda, mu, age=100, mrca=TRUE)[[1]];
## "(((((t35:2.336518061,t32:2.336518061):28.95257479,t10:31.28909285):8.654086516,t18:39.94317937):52.28906298,(((t47:31.00652286,t9:31.00652286):50.20634817,(t43:15.06939472,t38:15.06939472):66.14347631):10.61662549,(((((t20:20.94406932,t14:20.94406932):28.09437292,t19:49.03844224):31.16991698,(t24:54.88723469,(((t50:2.534803909,t8:2.534803909):35.18774941,t25:37.72255332):15.41876911,(t12:42.01655137,t5:42.01655137):11.12477106):1.745912255):25.32112453):2.610368667,t37:82.81872788):6.617999642,(t42:81.65977864,(t13:5.88018515,t41:5.88018515):75.77959349):7.776948892):2.392768999):0.4027458181):7.767757656,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0;"
tr <- paintSubTree(tr, node=78, state=2, stem=TRUE)
tr <- paintSubTree(tr, node=60, state=3, stem=TRUE)
plotSimmap(tr)

paste(tr$tip.label, collapse=",")

set.seed(1)
sigmas <- list(black=0.01, red=0.05, green=0.1)
dat <- mvSIM(tr, nsim=1, model="BMM", param=list(ntraits=1, sigma=sigmas, theta=0.0))
res <- mvBM(tr, dat, model="BMM")

res$sigma # 0.0082732 0.0509627 0.0953709
res$theta # 0.1513956
res$LogLik # -72.44127

## (5) 'BMMVNShiftLikelihoodOneTraitTest3'
