
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
### (1.1) Simulating one trait under Brownian motion (BM), writing .xmls from template and .sh scripts.

We use an exponential prior for the evolutionary rate (sigma^2), with rate 5. We use a normal prior for the mean of the process (=root value) with mean=0.0 and stdev=2.0.

```
cd /path/to/calibrated_validation/
mkdir BMMVNOneTrait_ultra_xmls/
mkdir BMMVNOneTrait_ultra_shellscripts/

Rscript BM_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVN '00:30:00' BMMVNLikelihoodOneTrait_fixedtree_template.xml ultra BMMVNLikelihoodOneTrait_fixedtree_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar
```

Here we simulate 100 data sets, on a tree with 50 species, set the job prefixes to 'BMMVN', and tell the cluster the jobs should take 30 minutes each. Here, "ultra" because I'm using an ultrametric tree, and the other arguments are self-explanatory. You can then simulate under a nonultrametric tree by just replacing the "tr" variable (and changing other arguments accordingly (see below). I hard-coded this in the script -- just uncomment it, or use whatever tree you want.

```
cd /path/to/calibrated_validation/
mkdir BMMVNOneTrait_nonultra_xmls/
mkdir BMMVNOneTrait_nonultra_shellscripts/

Rscript BM_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVN '00:30:00' BMMVNLikelihoodOneTrait_fixedtree_nonultra_template.xml nonultra BMMVNLikelihoodOneTrait_fixedtree_nonultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar
```

### (1.2) Plotting the mean posterior of the two BM parameters against the true values

```
Rscript BM_calibrated_postMCMC_plots.R ./ BMMVNLikelihoodOneTrait_fixedtree.RData 100 BMMVN ultra
Rscript BM_calibrated_postMCMC_plots.R ./ BMMVNLikelihoodOneTrait_fixedtree_nonultra.RData 100 BMMVN nonultra
```

## (2) Brownian motion with multiple rates (BMShift) multivariate normal likelihood 
### (2.1) Simulating one trait, but just one rate on whole tree under Brownian motion, writing .xmls from template and .sh scripts

We use the same priors as in 1.1, and do both ultrametric and nonultrametric trees.

```
cd /path/to/calibrated_validation/
mkdir BMMVNShiftOneRateOneTrait_ultra_xmls/
mkdir BMMVNShiftOneRateOneTrait_ultra_shellscripts/

Rscript BMShiftOneRate_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVNShiftOneRate '00:45:00' BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_ultra_template.xml ultra BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_ultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar

# remember to un/comment the tree specification lines in BMShiftOneRate_calibrated_simulation.R

Rscript BMShiftOneRate_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVNShiftOneRate '00:45:00' BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_nonultra_template.xml nonultra BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_nonultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar
```

### (2.2) Plotting the mean posterior of the two BM parameters against the true values

```
Rscript BMShiftOneRate_calibrated_postMCMC_plots.R ./ BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_ultra.RData 100 BMMVNShiftOneRate ultra
Rscript BMShiftOneRate_calibrated_postMCMC_plots.R ./ BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_nonultra.RData 100 BMMVNShiftOneRate nonultra
```

## (2.3) Simulating one trait, with three rates on whole tree under Brownian motion, writing .xmls from template and .sh scripts 

Note that BEAST2 will number internal nodes according to either ASCII-sorted names, or following the order of the provided TaxonSet. This means that the setting of different rates will only match that of BMMVNShiftLikelihoodOneTraitTest5 if we let ASCII-betical sorting happen (i.e., by not providing a TaxonSet). So template here is different from the usual.

Also, note that inside BMShiftThreeRates_calibrated_simulation.R, mvMORPH seems to call the third sigma^2 (the green color) "sigma[[2]]". So I had to adjust these accordingly both in this script and in the plotting script.

```
cd /path/to/calibrated_validation/
mkdir BMMVNShiftThreeRatesOneTrait_ultra_xmls/
mkdir BMMVNShiftThreeRatesOneTrait_ultra_shellscripts/
mkdir BMMVNShiftThreeRatesOneTrait_nonultra_xmls/
mkdir BMMVNShiftThreeRatesOneTrait_nonultra_shellscripts/

Rscript BMShiftThreeRates_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVNShiftThreeRates '01:30:00' BMMVNShiftLikelihoodThreeRatesOneTrait_fixedtree_ultra_template.xml ultra BMMVNShiftLikelihoodThreeRatesOneTrait_fixedtree_ultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar

# remember to un/comment the tree specification lines in BMShiftThreeRates_calibrated_simulation.R

Rscript BMShiftThreeRates_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVNShiftThreeRates '01:30:00' BMMVNShiftLikelihoodThreeRatesOneTrait_fixedtree_nonultra_template.xml nonultra BMMVNShiftLikelihoodThreeRatesOneTrait_fixedtree_nonultra_ /nesi/project/nesi00390/fkmendes/contraband/calibrated_validation /nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/contraband.jar
```

## (2.4) Plotting the mean posterior of the four BM parameters (three sigmas^2 and root value)

```
Rscript BMShiftThreeRates_calibrated_postMCMC_plots.R ./ BMMVNShiftLikelihoodThreeRatesOneTrait_fixedtree_ultra.RData 100 BMMVNShiftThreeRates ultra
Rscript BMShiftThreeRates_calibrated_postMCMC_plots.R ./ BMMVNShiftLikelihoodThreeRatesOneTrait_fixedtree_nonultra.RData 100 BMMVNShiftThreeRates nonultra
```

## (3) Brownian motion (BM) pruning likelihood
### (3.1) Simulating one trait under Brownian Motion (BM), writing .xmls from template and .sh scripts

We use the same tree and data sets from the BMMVN calibration.

```
cd /path/to/calibrated_validation/
mkdir BMPrune_ultra_xmls/
mkdir BMPrune_ultra_shellscripts/

Rscript BM_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMPrune '00:15:00' BMPruneLikelihoodOneTrait_fixedtree_template.xml

```

Then we uncomment the hard-coded non-ultrametric tree, and run simulations on it as well.

```
cd /path/to/calibrated_validation/
mkdir BMPrune_nonultra_xmls/
mkdir BMPrune_nonultra_shellscripts/

Rscript BM_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMPrune '00:15:00' BMPruneLikelihoodOneTrait_fixedtree_nonultra_template.xml
```

### (3.2) Plotting the mean posterior of the two BM parameters against the true values

```
Rscript BM_calibrated_postMCMC_plots.R ./ BMPruneLikelihoodOneTrait_fixedtree.RData 100 BMPrune ultra
Rscript BM_calibrated_postMCMC_plots.R ./ BMPruneLikelihoodOneTrait_fixedtree_nonultra.RData 100 BMPrune nonultra
```

## (4) Ornstein-Uhlenbeck (OU) multivariate normal likelihood
### (4.1) Simulating one trait under the vanilla (1 optimum for the whole tree) Ornstein-Uhlenbeck (OU), writing .xmls from template and .sh scripts.

We use an exponential prior (rate=5) for the evolutionary rate, sigma^2. We use a normal prior for the root value (mean=0.0 and stdev=2.0) and optimum (mean=1.0, stdev=2.0), and a lognormal prior for alpha (mean = stdev = 1.0).

```
cd /path/to/calibrated_validation/
mkdir OUVanillaMVNOneTrait_ultra_xmls/
mkdir OUVanillaMVNOneTrait_ultra_shellscripts/

Rscript OUVanilla_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 OUVanillaMVN '00:45:00' OUVanillaMVNLikelihoodOneTrait_fixedtree_template.xml ultra OUVanillaMVNLikelihoodOneTrait_fixedtree_ /nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/ /nesi/project/nesi00390/fkmendes/contraband/contraband.jar
```

We follow the same procedure as in BM. There is a hard-coded tree (you must uncomment it) that just makes some of the branch lengths very short, but you can edit that to be whatever tree you want. We use the same priors, but simulate along this non-ultrametric tree instead.

```
cd /path/to/calibrated_validation/
mkdir OUVanillaMVNOneTrait_nonultra_xmls/
mkdir OUVanillaMVNOneTrait_nonultra_shellscripts/

Rscript OUVanilla_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 OUVanillaMVN '00:45:00' OUVanillaMVNLikelihoodOneTrait_fixedtree_nonultra_template.xml nonultra OUVanillaMVNLikelihoodOneTrait_fixedtree_nonultra_ /nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/ /nesi/project/nesi00390/fkmendes/contraband/contraband.jar
```

### (4.2) Plotting the mean posterior of the four vanilla OU parameters against the true values

```
Rscript OUVanilla_calibrated_postMCMC_plots.R ./ OUVanillaMVNLikelihoodOneTrait_fixedtree.RData 100 OUVanillaMVN ultra
Rscript OUVanilla_calibrated_postMCMC_plots.R ./ OUVanillaMVNLikelihoodOneTrait_fixedtree_nonultra.RData 100 OUVanillaMVN nonultra
```
