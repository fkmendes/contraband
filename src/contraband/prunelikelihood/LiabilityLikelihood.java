package contraband.prunelikelihood;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import contraband.utils.PruneLikelihoodUtils;
import outercore.parameter.KeyRealParameter;

// This class implements likelihood of combined morphological data which has been transformed to continuous variables.
public class LiabilityLikelihood extends BMPruneLikelihood{
    final public Input<RealParameter> binaryLiabilitiesInput = new Input<>("binaryLiability", "Continuous random variables for  binary discrete traits.");
    final public Input<RealParameter> orderedLiabilitiesInput = new Input<>("orderedLiability", "Continuous random variables for ordered discrete traits.");

    private Tree tree;
    private int contTraitNr = 0;
    private int binaryTraitNr = 0;
    private int orderedTraitNr = 0;
    private int totalTraitNr;

    private int nSpecies;

    private KeyRealParameter contTraits;
    private double[] binaryLiabilities;
    private double[] orderedLiabilities;
    private double[] traitValuesArr;



    @Override
    public void initAndValidate() {
        // get the tree
        tree = treeInput.get();
        nSpecies = treeInput.get().getLeafNodeCount();

        // get liabilities for discrete traits
        if(binaryLiabilitiesInput.get() != null) {
            binaryTraitNr = binaryLiabilitiesInput.get().getDimension() / nSpecies;
            binaryLiabilities = binaryLiabilitiesInput.get().getDoubleValues();
        }
        if(orderedLiabilitiesInput.get() != null) {
            orderedTraitNr = orderedLiabilitiesInput.get().getDimension() / nSpecies;
            orderedLiabilities = orderedLiabilitiesInput.get().getDoubleValues();
        }

        // get continuous trait values
        if(traitsValuesInput.get() != null) {
            contTraits = traitsValuesInput.get();
            contTraitNr = contTraits.getMinorDimension1();
        }

        // get total number of traits
        totalTraitNr = binaryTraitNr + orderedTraitNr + contTraitNr;
        if (totalTraitNr == 0) {
            throw new RuntimeException("LiabilityLikelihood::At least one kind of data should be input.");
        }
        setNTraits(totalTraitNr);

        super.initAndValidate();
    }

    private void populateCombinedTraitValuesArr(double[] liabilities, double[] traitValuesArr, int index, int traitNr, int totalNr, int speciesNr){
        for (int i = 0; i < speciesNr; i++) {
            System.arraycopy(liabilities, i * traitNr, traitValuesArr, i * totalNr + index, traitNr);
        }
    }

    @Override
    protected void populateTraitData(){
        /*
         * populate combined trait values in an array
         */
        traitValuesArr = new double[nSpecies * totalTraitNr];
        int index = 0;
        // (1) populate continuous trait values
        if(traitsValuesInput.get() != null) {
            PruneLikelihoodUtils.populateTraitValuesArr(contTraits, tree, contTraitNr, traitValuesArr);
            index += contTraitNr;
        }
        // (2) populate binary discrete trait values
        if(binaryLiabilitiesInput.get() != null){
            populateCombinedTraitValuesArr(binaryLiabilities, traitValuesArr, index, binaryTraitNr, totalTraitNr, nSpecies);
            index += binaryTraitNr;
        }
        // (3) populate ordered discrete trait values
        if(orderedLiabilitiesInput.get() != null) {
            populateCombinedTraitValuesArr(orderedLiabilities, traitValuesArr, index, orderedTraitNr, totalTraitNr, nSpecies);
            index += orderedTraitNr;
        }
    }
}
