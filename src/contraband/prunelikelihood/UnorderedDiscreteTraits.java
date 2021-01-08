package contraband.prunelikelihood;

public class UnorderedDiscreteTraits extends ThresholdModel {

    private int nSpecies;
    private int nTraits;
    private int[] nrOfStates;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        nSpecies = getSpeciesNr();
        nTraits = getTraitNr();

        nrOfStates = new int[nTraits];
        populateTraitStates(nrOfStates, nTraits);


    }

    @Override
    public double calculateLogP() {
        return logP;
    }

    private int getLiabilityDimension(int[] states) {
        int sum = 0;
        for (int i : states){
            sum += (i-1);
        }
        return sum;
    }

}
