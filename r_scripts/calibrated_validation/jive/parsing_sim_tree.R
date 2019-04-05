library(ape)
library(phytools)

n.spp <- 100

# tree
tr <- read.tree("tree_n-100.tre")
taxon.strs.4.template <- paste(paste0("<taxon id=\"", tr$tip.label, "\" spec=\"Taxon\"/>"), collapse="\n                  ")

## getting species names for xml
cat(paste(tr$tip.label,collapse=","))

# trait data
traits <- read.table("traits_100sp_M-BM0.5_V-BM0.1.txt", header=TRUE)
names(traits)[3:ncol(traits)] <- paste0("sample",1:n.spp)
row.names(traits)[1:n.spp]

samples.4.xml <- rep("", n.spp)
for (i in 1:n.spp) {
    this.sp = row.names(traits)[i]
    this.sample = traits[i,][3:ncol(traits)]
    this.sample = this.sample[!is.na(this.sample)]
    samples.4.xml[i] = paste0(this.sp, "=", paste(this.sample, collapse=","))
}
cat(paste(samples.4.xml, collapse="|")) # for unit test
cat(paste(samples.4.xml, collapse="&#124;")) # for xml

## initial values for HyperMeanRootState and HyperMeanSigmaSq
hmrs <- runif(1, -1E4, 1E4)
hmss <- rgamma(1, 1.1, scale=5)

## initial values for JiveMeans
cat(paste(fastBM(tr, sig2=hmss, a=hmrs)))

## initial values for HyperVarMean and HyperVarSigmaSq
hvrs <- runif(1, -10, 10)
hvss <- rgamma(1, 1.1, scale=5)

## initial values for JiveLogVars
cat(paste(fastBM(tr, sig2=hvss, a=hvrs)))

## grabbing JiveMeans and JiveLogVars for unit test
cat(paste(traits[,1][1:100], collapse=",")) # logvars
cat(paste(traits[,2][1:100], collapse=",")) # means

#########################################

## need to do a very simple case first
hmrs <- 0.0
hmss <- 1.0
hvrs <- 0.5
hvss <- 0.1

tr <- read.tree(text="((sp1:1.0,sp2:1.0):2.0,sp3:3.0);")

hms <- fastBM(tr, sig2=hmss, a=hmrs)
hvs <- fastBM(tr, sig2=hvss, a=hvrs)
sp1 <- rnorm(3, mean=hms[1], sd=sqrt(hvs)[1])
sp2 <- rnorm(3, mean=hms[2], sd=sqrt(hvs)[2])
sp3 <- rnorm(3, mean=hms[3], sd=sqrt(hvs)[3])

paste0("sp1=",paste(sp1, collapse=","))
paste0("sp2=",paste(sp2, collapse=","))
paste0("sp3=",paste(sp3, collapse=","))
