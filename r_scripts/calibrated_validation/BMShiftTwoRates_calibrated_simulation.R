library(TreeSim)
library(phytools)
library(stringr)
source("calibrated_validation_utils.R")

args = commandArgs(trailingOnly=TRUE)

### SCRIPT FLAGS AND PATH VARIABLES ###

simulate <- TRUE
## simulate <- args[1]
write.xmls <- TRUE
## write.xmls <- args[2]
write.shellscripts <- TRUE
## write.shellscripts <- args[3]
cal.validation.folder <- "./"
## cal.validation.folder <- args[4]
n.sim <- 100
## n.sim <- as.numeric(args[5])
n.spp <- 50
## n.spp <- as.numeric(args[6])
job.prefix <- "BMMVNShiftTwoRates"
## job.prefix <- args[7] # e.g., "BMMVNShiftTwoRates" or "BMPruneShiftTwoRates"
time.needed <- "01:30:00"
## time.needed <- args[8] # e.g., "01:30:00"
template.name <- "BMMVNShiftLikelihoodTwoRatesOneTrait_fixedtree_ultra_template.xml"
## template.name <- args[9]
tree.type <- "ultra"
## tree.type <- args[10]
xmlfolder.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_xmls/")
xml.file.prefix <- "BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_ultra_"
## xml.file.prefix <- args[11]
shell.scripts.path <- paste0(cal.validation.folder, job.prefix, "OneTrait_", tree.type, "_shellscripts/")
template.path <- paste0(cal.validation.folder, template.name)
rdata.path <- gsub("_template.xml", ".RData", template.path)

# cluster stuff
cluster.validation.folder <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation"
## cluster.validation.folder <- args[12]
xml.file.path <- paste0(cluster.validation.folder, "BMMVNShiftTwoRatesOneTrait_", tree.type, "_xmls/")
res.path <- paste0(cluster.validation.folder, "BMMVNShiftTwoRatesOneTrait_", tree.type, "_results/")
jar.path <- "/nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/contraband.jar"
## jar.path <- args[13]

n.param <- 3
sigma.rate <- 5 # exponential prior
x0.mean <- 0.0 # normal prior
x0.sd <- 2.0

############# DOING STUFF #############

## in ape, branch numbering is done in postorder order: to see the order of branches, you can do tr <- reorder(tr, "po"); plor(tr); edgelabels()

set.seed(123)

## simulating birth-death tree
lambda <- rexp(1, rate=80) # lambda from exponential (mean = 1/80)
mu <- rexp(1, rate=100) # mu from exponential (mean = 1/100)
tr <- sim.bd.taxa.age(n.spp, 1, lambda, mu, age=100, mrca=TRUE)[[1]];
## "(((((t35:2.336518061,t32:2.336518061):28.95257479,t10:31.28909285):8.654086516,t18:39.94317937):52.28906298,(((t47:31.00652286,t9:31.00652286):50.20634817,(t43:15.06939472,t38:15.06939472):66.14347631):10.61662549,(((((t20:20.94406932,t14:20.94406932):28.09437292,t19:49.03844224):31.16991698,(t24:54.88723469,(((t50:2.534803909,t8:2.534803909):35.18774941,t25:37.72255332):15.41876911,(t12:42.01655137,t5:42.01655137):11.12477106):1.745912255):25.32112453):2.610368667,t37:82.81872788):6.617999642,(t42:81.65977864,(t13:5.88018515,t41:5.88018515):75.77959349):7.776948892):2.392768999):0.4027458181):7.767757656,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0;"

## tr <- read.tree(text="(((((t35:0.1,t32:0.1):0.1,t10:0.1):0.1,t18:0.1):0.1,(((t47:0.1,t9:0.1):0.1,(t43:0.1,t38:0.1):0.1):0.1,(((((t20:0.1,t14:0.1):0.1,t19:0.1):0.1,(t24:0.1,(((t50:0.1,t8:0.1):0.1,t25:0.1):0.1,(t12:0.1,t5:0.1):0.1):0.1):0.1):0.1,t37:0.1):0.1,(t42:0.1,(t13:0.1,t41:0.1):0.1):0.1):0.1):0.1):0.1,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0;")

tr <- paintSubTree(tr, node=78, state=2, stem=TRUE)
tr <- paintSubTree(tr, node=60, state=3, stem=TRUE)
## plotSimmap(tr)
rate.assignments <- # gotta run JUnit test in contraband to see how to assign rates

## simulating quant trait data sets
sigmas1 <- rexp(n.sim, rate=sigma.rate); # sigmas1[1] # 0.02914
sigmas2 <- rexp(n.sim, rate=sigma.rate); # sigmas2[1] # 0.06309
sigmas3 <- rexp(n.sim, rate=sigma.rate); # sigmas3[1] # 0.34130
x0s <- rnorm(n.sim, mean=x0.mean, sd=x0.sd); # x0s[1] # 2.44268
datasets <- vector("list", n.sim) # storing sims
mles <- data.frame(matrix(NA,100,n.param))

