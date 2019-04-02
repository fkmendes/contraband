# contraband: Calibrated validation of BM-related classes

This file will guide you through reproducing the simulations and graphs for conducting calibrated validation of the relevant likelihoods in the contraband package.
You are going to need the following R packages:

* mvMORPH
* TreeSim
* phytools
* stringr
* ggplot2
* gtools
* sjPlot

Also, make sure the file calibrated_validation_utils.R is inside the calibrated_validation/ folder, or whatever folder you use for holding all files.

## (1) Brownian motion (BM) multivariate normal likelihood
### (1.1) Simulating one trait under Brownian Motion (BM), writing .xmls from template and .sh scripts.

We use an exponential prior for the evolutionary rate (sigma^2), with rate 5. We use a normal prior for the mean of the process (=root value) with mean=0.0 and stdev=2.0.

```
cd /path/to/calibrated_validation/
mkdir BMMVNOneTrait_ultra_xmls/
mkdir BMMVNOneTrait_ultra_shellscripts/

Rscript BM_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVN '00:30:00' BMMVNLikelihoodOneTrait_fixedtree_template.xml ultra BMMVNLikelihoodOneTrait_fixedtree_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar
```

Here we simulate 100 data sets, on a tree with 50 species, set the job prefixes to 'BMMVN', and tell the cluster the jobs should take 30 minutes each. Here, "ultra" because I'm using an ultrametric tree, and the other arguments are self-explanatory. You can then simulate under a nonultrametric tree by just replacing the "tr" variable (and changing other arguments accordingly (see below). I hard-coded this in the script -- just uncomment it, or use whatever tree you want.

```
mkdir BMMVNOneTrait_nonultra_xmls/
mkdir BMMVNOneTrait_nonultra_shellscripts/

Rscript BM_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVN '00:30:00' BMMVNLikelihoodOneTrait_fixedtree_nonultra_template.xml nonultra BMMVNLikelihoodOneTrait_fixedtree_nonultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar
```

### (1.2) Plotting the mean posterior of the two BM parameters against the true values

```
Rscript BM_calibrated_postMCMC_plots.R ./ BMMVNLikelihoodOneTrait_fixedtree_template.xml BMMVN ultrametric
Rscript BM_calibrated_postMCMC_plots.R ./ BMMVNLikelihoodOneTrait_fixedtree_nonultra_template.xml BMMVN nonultrametric
```

## (2) Brownian motion (BM) pruning likelihood
### (2.1) Simulating one trait under Brownian Motion (BM), writing .xmls from template and .sh scripts.

We use the same tree and data sets from the BMMVN calibration.

```
cd /path/to/calibrated_validation/
mkdir BMPrune_xmls/
mkdir BMPrune_shellscripts/

Rscript BM_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMPrune '00:10:00' BMPruneLikelihoodOneTrait_fixedtree_template.xml
```

### (2.2) Plotting the mean posterior of the two BM parameters against the true values

```
Rscript BM_calibrated_postMCMC_plots.R ./ BMPruneLikelihoodOneTrait_fixedtree_template.xml BMPrune
```