package contraband.prunelikelihood;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import outercore.parameter.KeyRealParameter;

// This class implements likelihood of combined morphological data which has been transformed to continuous variables.
public class LiabilityLikelihood extends BMPruneLikelihood{
    final public Input<RealParameter> binaryLiabilitiesInput = new Input<>("binaryLiability", "Continuous random variables for  binary discrete traits.");
    final public Input<RealParameter> orderedLiabilitiesInput = new Input<>("orderedLiability", "Continuous random variables for ordered discrete traits.");
    final public Input<UnorderedDiscreteTraits> unorderedDiscreteTraitsInput = new Input<>("unorderedDiscreteTraits", "Object for unordered discrete traits.");

    private Tree tree;
    private int contTraitNr = 0;
    private int binaryTraitNr = 0;
    private int orderedTraitNr = 0;
    private int unorderedLiabilityNr = 0;
    private int totalTraitNr;

    private int nSpecies;
    private KeyRealParameter contTraits;
    private double[] traitValuesArr;

    private int binaryTraitIndex;
    private int orderedTraitIndex;
    private  int unorderedLiabilityIndex;

    @Override
    public void initAndValidate() {
        // get the tree
        tree = treeInput.get();
        nSpecies = treeInput.get().getLeafNodeCount();

        // get liabilities for discrete traits
        if(binaryLiabilitiesInput.get() != null) {
            binaryTraitNr = binaryLiabilitiesInput.get().getDimension() / nSpecies;
        }
        if(orderedLiabilitiesInput.get() != null) {
            orderedTraitNr = orderedLiabilitiesInput.get().getDimension() / nSpecies;
        }
        if(unorderedDiscreteTraitsInput.get() != null){
            UnorderedDiscreteTraits unorderedDiscreteTraits = unorderedDiscreteTraitsInput.get();
            unorderedLiabilityNr = unorderedDiscreteTraits.getNrOfLiabilities();
        }

        // get continuous trait values
        if(traitsValuesInput.get() != null) {
            contTraits = traitsValuesInput.get();
            contTraitNr = contTraits.getMinorDimension1();
        }

        // get total number of traits
        totalTraitNr = binaryTraitNr + orderedTraitNr + unorderedLiabilityNr + contTraitNr;
        if (totalTraitNr == 0) {
            throw new RuntimeException("LiabilityLikelihood::At least one kind of data should be input.");
        }
        setNTraits(totalTraitNr);

        super.initAndValidate();

        setTraitValuesArr(traitValuesArr);
    }

    // add liabilities of discrete traits to the traitValuesArr
    private void populateCombinedTraitValuesArr(double[] liabilities, double[] traitValuesArr, int index, int traitNr, int totalNr, int speciesNr){
        for (int i = 0; i < speciesNr; i++) {
            System.arraycopy(liabilities, i * traitNr, traitValuesArr, i * totalNr + index, traitNr);
        }
    }

    // add continuous trait values to the traitValuesArr
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
        binaryTraitIndex = index;

        // (2) populate binary discrete trait values
        if(binaryLiabilitiesInput.get() != null){
            populateCombinedTraitValuesArr(binaryLiabilitiesInput.get().getDoubleValues(), traitValuesArr, index, binaryTraitNr, totalTraitNr, nSpecies);
            index += binaryTraitNr;
        }
        orderedTraitIndex = index;

        // (3) populate ordered discrete trait values
        if(orderedLiabilitiesInput.get() != null) {
            populateCombinedTraitValuesArr(orderedLiabilitiesInput.get().getDoubleValues(), traitValuesArr, index, orderedTraitNr, totalTraitNr, nSpecies);
            index += orderedTraitNr;
        }
        unorderedLiabilityIndex = index;

        // (4) populate unordered discrete trait values
        if(unorderedDiscreteTraitsInput.get() != null){
            double[] unorderedLiabilities = unorderedDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(unorderedLiabilities, traitValuesArr, index, unorderedLiabilityNr, totalTraitNr, nSpecies);
        }

    }

    public double[] getCombinedTraitDataArr() { return traitValuesArr; }

    @Override
    public double calculateLogP() {
        updateTraitValuesArr();

        logP = super.calculateLogP();

        return logP;
    }

    private void updateTraitValuesArr(){
        // (1) update binary discrete trait values
        if(binaryLiabilitiesInput.get() != null){
            populateCombinedTraitValuesArr(binaryLiabilitiesInput.get().getDoubleValues(), traitValuesArr, binaryTraitIndex, binaryTraitNr, totalTraitNr, nSpecies);
        }

        // (2) update ordered discrete trait values
        if(orderedLiabilitiesInput.get() != null) {
            populateCombinedTraitValuesArr(orderedLiabilitiesInput.get().getDoubleValues(), traitValuesArr, orderedTraitIndex, orderedTraitNr, totalTraitNr, nSpecies);
        }

        // (3) update unordered discrete trait values
        if(unorderedDiscreteTraitsInput.get() != null){
            double[] unorderedLiabilities = unorderedDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(unorderedLiabilities, traitValuesArr, unorderedLiabilityIndex, unorderedLiabilityNr, totalTraitNr, nSpecies);
        }

        setTraitValuesArr(traitValuesArr);
    }
}
