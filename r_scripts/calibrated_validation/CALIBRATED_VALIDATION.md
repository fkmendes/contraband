
# contraband: Calibrated validation of BM-related classes

This file will guide you through reproducing the simulations and graph-plotting we did for our well-calibrated validation. Note that first we fix the tree to verify the BM-related likelihoods by themselves, and but we also validate with respect to a tree prior.

You are going to need the following R packages:

* mvMORPH
* TreeSim
* phytools
* stringr
* ggplot2
* gtools
* sjPlot

Also, make sure the file calibrated_validation_utils.R is inside the calibrated_validation/ folder, or whatever folder you use for holding all files.

---

## Fixed tree calibrated validation
## (1) Brownian motion (BM) multivariate normal likelihood
### (1.1) Simulating one trait under Brownian motion (BMMVN class), writing .xmls from template and .sh scripts.

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

## (2) Brownian motion with potentially multiple rates (BMMVNShift class) multivariate normal likelihood 
### (2.1) Simulating one trait, but just one rate on whole tree under Brownian motion, writing .xmls from template and .sh scripts

We use the same priors as in 1.1, and do ultrametric and nonultrametric trees like before.

```
cd /path/to/calibrated_validation/
mkdir BMMVNShiftOneRateOneTrait_ultra_xmls/
mkdir BMMVNShiftOneRateOneTrait_ultra_shellscripts/
mkdir BMMVNShiftOneRateOneTrait_nonultra_xmls/
mkdir BMMVNShiftOneRateOneTrait_nonultra_shellscripts/

Rscript BMShiftOneRate_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVNShiftOneRate '00:45:00' BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_ultra_template.xml ultra BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_ultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar

# remember to un/comment the tree specification lines in BMShiftOneRate_calibrated_simulation.R

Rscript BMShiftOneRate_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVNShiftOneRate '00:45:00' BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_nonultra_template.xml nonultra BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_nonultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar
```

### (2.2) Plotting the mean posterior of the two BM parameters against the true values

```
Rscript BMShiftOneRate_calibrated_postMCMC_plots.R ./ BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_ultra.RData 100 BMMVNShiftOneRate ultra
Rscript BMShiftOneRate_calibrated_postMCMC_plots.R ./ BMMVNShiftLikelihoodOneRateOneTrait_fixedtree_nonultra.RData 100 BMMVNShiftOneRate nonultra
```

### (2.3) Simulating one trait, with three rates on whole tree under Brownian motion, writing .xmls from template and .sh scripts 

Note that BEAST2 will number internal nodes according to either ASCII-sorted names, or following the order of the provided TaxonSet. This means that the setting of different rates will only match that of BMMVNShiftLikelihoodOneTraitTest5 if we let ASCII-betical sorting happen (i.e., by not providing a TaxonSet). So template here is different from the usual.

Also, note that inside BMShiftThreeRates_calibrated_simulation.R, mvMORPH seems to call the third sigma^2 (the green color) "sigma[[2]]". So I had to adjust these accordingly both in this script and in the plotting script.

As before, we validate on ultrametric and nonultrametric trees.

```
cd /path/to/calibrated_validation/
mkdir BMMVNShiftThreeRatesOneTrait_ultra_xmls/
mkdir BMMVNShiftThreeRatesOneTrait_ultra_shellscripts/
mkdir BMMVNShiftThreeRatesOneTrait_nonultra_xmls/
mkdir BMMVNShiftThreeRatesOneTrait_nonultra_shellscripts/

Rscript BMShiftThreeRates_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVNShiftThreeRates '01:00:00' BMMVNShiftLikelihoodThreeRatesOneTrait_fixedtree_ultra_template.xml ultra BMMVNShiftLikelihoodThreeRatesOneTrait_fixedtree_ultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar 

# remember to un/comment the tree specification lines in BMShiftOneRate_calibrated_simulation.R

Rscript BMShiftThreeRates_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVNShiftThreeRates '01:00:00' BMMVNShiftLikelihoodThreeRatesOneTrait_fixedtree_nonultra_template.xml nonultra BMMVNShiftLikelihoodThreeRatesOneTrait_fixedtree_nonultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar 
```

### (2.4) Plotting the mean posterior of the four BM parameters (three sigmas^2 and root value)

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
mkdir OUMVNVanillaOneTrait_ultra_xmls/
mkdir OUMVNVanillaOneTrait_ultra_shellscripts/
mkdir OUMVNVanillaOneTrait_nonultra_xmls/
mkdir OUMVNVanillaOneTrait_nonultra_shellscripts/

Rscript OUVanilla_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 OUMVNVanilla '01:15:00' OUMVNLikelihoodVanillaOneTrait_fixedtree_ultra_template.xml ultra OUMVNLikelihoodVanillaOneTrait_fixedtree_ultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar

# remember to un/comment the tree specification lines in OUVanilla_calibrated_simulation

Rscript OUVanilla_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 OUMVNVanilla '01:20:00' OUMVNLikelihoodVanillaOneTrait_fixedtree_nonultra_template.xml nonultra OUMVNLikelihoodVanillaOneTrait_fixedtree_nonultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar
```

### (4.2) Plotting the mean posterior of the four vanilla OU parameters against the true values

```
Rscript OUVanilla_calibrated_postMCMC_plots.R ./ OUMVNLikelihoodVanillaOneTrait_fixedtree_ultra.RData 100 OUMVNVanilla ultra
Rscript OUVanilla_calibrated_postMCMC_plots.R ./ OUMVNLikelihoodVanillaOneTrait_fixedtree_nonultra.RData 100 OUMVNVanilla nonultra
```

### (4.3) Simulating one trait with three adaptive optima under Orstein-Uhlenbeck (OU), writing .xmls from template and .sh scripts.

The priors we use are the same as those in (4.1).

```
cd /path/to/calibrated_validation/
mkdir OUMVNThreeOptOneTrait_ultra_xmls/
mkdir OUMVNThreeOptOneTrait_ultra_shellscripts/
mkdir OUMVNThreeOptOneTrait_nonultra_xmls/
mkdir OUMVNThreeOptOneTrait_nonultra_shellscripts/

Rscript OUThreeOpt_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 OUMVNThreeOpt '01:45:00' OUMVNLikelihoodThreeOptOneTrait_fixedtree_ultra_template.xml ultra OUMVNThreeOptOneTrait_fixedtree_ultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar

# remember to un/comment the tree specification lines in OUThreeOpt_calibrated_simulation

Rscript OUThreeOpt_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 OUMVNThreeOpt '01:45:00' OUMVNLikelihoodThreeOptOneTrait_fixedtree_nonultra_template.xml nonultra OUMVNThreeOptOneTrait_fixedtree_nonultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar
```

### (4.4) Plotting the mean posterior of the six OU parameters (three optima) against the true values

```
Rscript OUThreeOpt_calibrated_postMCMC_plots.R ./ OUMVNLikelihoodThreeOptOneTrait_fixedtree_ultra.RData 100 OUMVNThreeOpt ultra 6 sigmasq,rv,theta1,theta2,theta3,alpha OUSigmaSq,OURootValue,OUTheta1,OUTheta3,OUTheta2,OUAlpha "expression(sigma^2),expression(y[0]),expression(theta[1]),expression(theta[2]),expression(theta[3]),expression(alpha)" "1/5,0.0,0.0,1.0" sigmasq.mle,rv.mle,theta1.mle,theta2.mle,theta3.mle,alpha.mle

Rscript OUThreeOpt_calibrated_postMCMC_plots.R ./ OUMVNLikelihoodThreeOptOneTrait_fixedtree_nonultra.RData 100 OUMVNThreeOpt nonultra 6 sigmasq,rv,theta1,theta2,theta3,alpha OUSigmaSq,OURootValue,OUTheta1,OUTheta3,OUTheta2,OUAlpha "expression(sigma^2),expression(y[0]),expression(theta[1]),expression(theta[2]),expression(theta[3]),expression(alpha)" "1/5,0.0,0.0,1.0" sigmasq.mle,rv.mle,theta1.mle,theta2.mle,theta3.mle,alpha.mle
```

---

## Sampling tree validation
## (5) Brownian motion with multiple rates (BMMVNShift class) multivariate normal likelihood 
## (5.1) Simulating one trait, but just one rate on whole tree under Brownian motion, writing .xmls from template and .sh scripts

We use the same priors as in 1.1, and do nonultrametric trees only using an FBD prior. The prior on FBD parameters are exponentials with rates 80, 100 and 150 (for lambda, mu and psi, respectively).

```
Rscript BMShiftOneRateFBD_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVNShiftOneRateFBD '03:30:00' BMMVNShiftLikelihoodOneRateFBDOneTrait_nonultra_template.xml nonultra "BMMVNShiftLikelihoodOneRateFBDOneTrait_nonultra_ /cluster/validation/folder/ /cluster/jarfilepath/contraband.jar
```

## (5.2) Plotting the mean posterior of the two BM parameters (rate and mean/ancestral value) against true values

```
Rscript BMShiftOneRateFBD_calibrated_postMCMC_plots.R ./ BMMVNShiftLikelihoodOneRateFBDOneTrait_nonultra.RData 100 BMMVNShiftOneRateFBD 2 sigmasq,mu rateValues,BMMean "expression(sigma^2),expression(y[0])" "1/5,0.0" sigmasq.mle,mu.mle
```

### (5.3) Simulating one trait, with two rates on whole tree under Brownian motion, writing .xmls from template and .sh scripts

The same simulations as in 5.1, but with two rates. We pick a random internal node and all its descendants to have the second rate.

```
Rscript BMShiftTwoRatesFBD_calibrated_simulation.R TRUE TRUE TRUE ./ 100 50 BMMVNShiftTwoRatesFBD '06:00:00' BMMVNShiftLikelihoodTwoRatesFBDOneTrait_nonultra_template.xml nonultra "BMMVNShiftLikelihoodTwoRatesFBDOneTrait_nonultra_ /nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/ /nesi/project/nesi00390/fkmendes/contraband/calibrated_validation/contraband.jar
``` 