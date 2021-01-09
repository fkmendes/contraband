package contraband.prunelikelihood;


import java.util.Arrays;

public class UnorderedDiscreteTraits extends ThresholdModel {

    private int nSpecies;
    private int nTraits;
    private int[] nrOfStates;
    private int nrOfLiabilities;
    private int[] liabilityIndex;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        nSpecies = getSpeciesNr();
        nTraits = getTraitNr();

        nrOfStates = new int[nTraits];
        populateTraitStates(nrOfStates, nTraits);
        nrOfLiabilities = getParameterDimension(nrOfStates);
        liabilitiesInput.get().setDimension(nrOfLiabilities * nSpecies);
        liabilityIndex = getParameterIndex(nrOfStates, nTraits);

        initiateLiabilities();
    }

    @Override
    protected void initiateLiabilities(){
       for(int traitIdx = 0; traitIdx < nTraits; traitIdx++){
                int liabilityNr = nrOfStates[traitIdx];
                for(int speciesIdx = 0; speciesIdx < nSpecies; speciesIdx ++){
                    int traitData = getTraitDataForSpecies(speciesIdx, traitIdx);
                    if(traitData == 0) {
                        // 0 = sup{0, liabilities}
                        // -> all liabilities = -1
                        for (int liabilityIdx = 0; liabilityIdx < liabilityNr; liabilityIdx++) {
                            liabilitiesInput.get().setValue(nrOfLiabilities * speciesIdx + liabilityIndex[traitIdx] +liabilityIdx, -1.0);
                        }
                    } else{
                        // the i-th liability = sup{0, liabilities}
                        // -> i-th liability = traitData
                        // -> other liabilities = -1
                        for (int liabilityIdx = 0; liabilityIdx < liabilityNr; liabilityIdx++) {
                            if(liabilityIdx == (traitData - 1)){
                                liabilitiesInput.get().setValue(nrOfLiabilities * speciesIdx + liabilityIndex[traitIdx] +liabilityIdx, traitData + 0.0);
                            } else {
                                liabilitiesInput.get().setValue(nrOfLiabilities * speciesIdx + liabilityIndex[traitIdx] +liabilityIdx, - 1.0);
                            }
                        }
                    }
            }
        }

    }

    @Override
    public double calculateLogP() {
        double[] liabilities = liabilitiesInput.get().getDoubleValues();
        for (int speciesIdx = 0; speciesIdx < nSpecies; speciesIdx ++){
            for(int traitIdx = 0; traitIdx < nTraits; traitIdx ++){
                int traitData = getTraitDataForSpecies(speciesIdx, traitIdx);
                // augmentedLiabilities = {0, liabilities}
                double[] augmentedLiabilities = new double[nrOfStates[traitIdx] + 1];
                System.arraycopy(liabilities, nrOfLiabilities * speciesIdx + liabilityIndex[traitIdx], augmentedLiabilities,1, nrOfStates[traitIdx]);
                // max = sup{augmentedLiabilities}
                double max = Arrays.stream(augmentedLiabilities).max().getAsDouble();
                // traitData-th element in {augmentedLiabilities} == max
                if(augmentedLiabilities[traitData] == max){
                    logP = 0.0;
                } else{
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }

        return logP;
    }

    public int getNrOfLiabilities(){ return nrOfLiabilities; }


}
