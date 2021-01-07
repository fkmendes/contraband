package contraband.prunelikelihood;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
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

        setTraitValuesArr(traitValuesArr);
    }

    private void populateCombinedTraitValuesArr(double[] liabilities, double[] traitValuesArr, int index, int traitNr, int totalNr, int speciesNr){
        for (int i = 0; i < speciesNr; i++) {
            System.arraycopy(liabilities, i * traitNr, traitValuesArr, i * totalNr + index, traitNr);
        }
    }

    private void populateContTraitArr(KeyRealParameter traitValues, Tree tree, int contTraitNr, int totalTraitNr, double[] traitValuesArr) {
        // according to node number of tips
        for (int i = 0; i < tree.getLeafNodeCount(); i ++) {
            // get all traits values for this species
            Double[] traitForSpecies = traitValues.getRowValues(tree.getNode(i).getID());
            for (int j= 0; j < contTraitNr; j ++) {
                // populate the traits one by one in an array
                traitValuesArr[i*totalTraitNr + j] = traitForSpecies[j];
            }
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
            populateContTraitArr(contTraits, tree, contTraitNr, totalTraitNr, traitValuesArr);
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

    public double[] getCombinedTraitDataArr() { return traitValuesArr; }
}
