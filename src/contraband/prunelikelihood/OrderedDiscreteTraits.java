package contraband.prunelikelihood;

import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.util.Arrays;


public class OrderedDiscreteTraits extends ThresholdModel {
    final public Input<RealParameter> thresholdsInput = new Input<>("threshold", "Thresholds for mapping ordered discrete traits.", Input.Validate.REQUIRED);

    private int nSpecies;
    private int nTraits;
    private int[] nrOfStates;
    private int[] thresholdIndex;

    @Override
    public void initAndValidate() {

        super.initAndValidate();
        nSpecies = getSpeciesNr();
        nTraits = getLiabilityNr();
        if (liabilitiesInput.get().getDimension() != nSpecies * nTraits){
            liabilitiesInput.get().setDimension(nSpecies * nTraits);
        }

        nrOfStates = new int[nTraits];
        populateTraitStates(nrOfStates, nTraits);

        initiateLiabilities();

        thresholdsInput.get().setDimension(getParameterDimension(nrOfStates));
        thresholdIndex = getParameterIndex(nrOfStates, nTraits);
        initiateThresholds();
    }

    @Override
    protected void initiateLiabilities(){
        /*
         * For example, if the data is [2 0 1],
         * the corresponding liabilities will be [3, -1, 1]
         */
        for(int i = 0; i < nSpecies; i++){
            for (int j = 0; j < nTraits; j ++){
                int traitValue = getTraitDataForSpecies(i, j);
                if(traitValue == 0) {
                    // trait value - 1
                    liabilitiesInput.get().setValue(nTraits * i + j, -1.0);
                }
                else if(traitValue == nrOfStates[j]){
                    // trait value + 1
                    liabilitiesInput.get().setValue(nTraits * i + j, traitValue + 1.0);
                }
                else {
                    // equivalent to its trait value
                    liabilitiesInput.get().setValue(nTraits * i + j, traitValue + 0.0);
                }
            }
        }
    }


    @Override
    public double calculateLogP() {
        double[] liabilities = liabilitiesInput.get().getDoubleValues();
        double[] thresholds = thresholdsInput.get().getDoubleValues();
        for (int traitIdx = 0; traitIdx < nTraits; traitIdx++){
            // sort the thresholds for this trait in ascending order
            double[] sortedThresholds = getThresholdsForTrait(traitIdx, thresholds);
            Arrays.sort(sortedThresholds);

            for (int speciesIdx = 0; speciesIdx < nSpecies; speciesIdx++){
                int traitData = getTraitDataForSpecies(speciesIdx, traitIdx);
                double liability = liabilities[speciesIdx * nTraits + traitIdx];
                boolean valid;
                if(traitData == 0){
                    // the first state must be smaller than the smallest threshold
                    valid = liability < sortedThresholds[0];
                }
                else if(traitData == nrOfStates[traitIdx]){
                    // the last state must be greater than the largest threshold
                    valid = liability > sortedThresholds[traitData-1];
                } else {
                    // the middle state mush be within the two thresholds
                    valid = (liability > sortedThresholds[traitData-1]) && (liability < sortedThresholds[traitData]);
                }

                // return likelihood
                if(valid){
                    logP = 0.0;
                } else {
                    logP = Double.NEGATIVE_INFINITY;
                }
            }
        }



        return logP;
    }

    private void initiateThresholds() {
        int idx = 0;
        for (int i = 0; i < nTraits; i++){
            int stateNr = nrOfStates[i];
            for (int j = 0; j < stateNr; j ++){
                double randomValue =  j + 0.5;
                thresholdsInput.get().setValue(idx, randomValue);
                idx ++;
            }
        }
    }

    public double[] getThresholds() {
        return thresholdsInput.get().getDoubleValues();
    }

    private double[] getThresholdsForTrait (int traitIdx, double[] thresholds){
        int nr = nrOfStates[traitIdx];
        double[] res = new double[nr];
        if (nr >= 0) System.arraycopy(thresholds, thresholdIndex[traitIdx], res, 0, nr);
        return res;
    }


}
