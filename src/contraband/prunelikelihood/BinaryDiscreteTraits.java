package contraband.prunelikelihood;


public class BinaryDiscreteTraits extends ThresholdModel{


    private int nSpecies;
    private int nTraits;

    @Override
    public void initAndValidate() {

        super.initAndValidate();

        nSpecies = getSpeciesNr();
        nTraits = getTraitNr();
        if (liabilitiesInput.get().getDimension() != nSpecies * nTraits){
            liabilitiesInput.get().setDimension(nSpecies * nTraits);
        }


    }

    @Override
    public double calculateLogP() {

        for (int i = 0; i < nSpecies; i++){
            for (int j = 0; j < nTraits; j ++){
                int traitValue = getDateForSpecies(i, j);
                double liability = liabilitiesInput.get().getArrayValue(nTraits * i + j);
                if (traitValue == 1 && liability > 0.0){
                    //  present trait data and positive liability value
                    logP = 0.0;
                } else if (traitValue == 0 && liability < 0.0){
                    //  absent trait data and negative liability value
                    logP = 0.0;
                } else {
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }
        return logP;
    }


}
