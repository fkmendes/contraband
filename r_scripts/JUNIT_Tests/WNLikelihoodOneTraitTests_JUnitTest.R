# author: Fabio K. Mendes
# This R script gives us the expected values for JUnit tests
# (1) 'WNLikelihoodOneTraitTest'

library(geiger)

## (1) 'WNLikelihoodOneTraitTest

tr <- read.tree(text="((sp1:1.0,sp2:1.0):1.0,sp3:2.0);")

nspp <- 3

set.seed(1)
dat <- data.frame(c(rnorm(nspp, mean=0, sd=1)))
rownames(dat) <- c("sp1", "sp2", "sp3")
res <- fitContinuous(tr, dat, model="white")
res$opt$sigsq # 0.9733856
res$opt$z0 # 0.3681067
res$opt$lnL # -4.216353

sum(log(dnorm(dat[,1], mean=0.3681067, sd=sqrt(0.9733856)))) ## -4.216453 (another way to do it...)
