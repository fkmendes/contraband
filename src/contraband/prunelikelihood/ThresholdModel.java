package contraband.prunelikelihood;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.DataType;
import beast.evolution.tree.Tree;
import java.util.List;
import java.util.Random;

// This class implements mappings from discrete traits to continuous random variables.
public abstract class ThresholdModel extends Distribution {
    final public Input<RealParameter> liabilitiesInput = new Input<>("liability", "Continuous random variables underneath observed traits.", Input.Validate.REQUIRED);
    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Input.Validate.REQUIRED);

    private int[] discreteDataArr;
    private List<Sequence> sequenceList;
    private List<String> taxaNames;
    private int nrOfSpecies;
    private DataType dataType;
    private int nrOfTraits;

    @Override
    public void initAndValidate() {

        Alignment data = dataInput.get();
        sequenceList = data.sequenceInput.get();

        taxaNames = data.getTaxaNames();
        nrOfTraits = data.getSiteCount();
        dataType = data.getDataType();

        Tree tree = treeInput.get();
        nrOfSpecies = tree.getLeafNodeCount();

        // populate an array of discrete data set
        discreteDataArr = new int[nrOfTraits * nrOfSpecies];
        populateDiscreteDataArray(tree, discreteDataArr);
    }

    // this method return an array of discrete trait values for species
    public int[] getAllDataForSpecies(String species){
        for (int i = 0; i < nrOfSpecies; i++) {
            if (taxaNames.get(i).equals(species)) {
                List<Integer> values = sequenceList.get(i).getSequence(dataType);
                int[] valuesArr = new int[nrOfTraits];
                for(int j = 0; j < nrOfTraits; j ++) {
                    valuesArr[j] = values.get(j).intValue();
                }
                return valuesArr;
            }
        }
        return null;
    }

    // this method return an array of traitIndex-th discrete trait values for all species
    public int[] getAllTraitData (int traitIndex) {
        int[] iThTraitValues = new int[nrOfSpecies];
        for (int i = 0; i < nrOfSpecies; i++) {
            List<Integer> values = sequenceList.get(i).getSequence(dataType);
            iThTraitValues[i] = values.get(traitIndex).intValue();
        }
        return iThTraitValues;
    }

    private void populateDiscreteDataArray (Tree tree, int[] arr){
        for (int i = 0; i < nrOfSpecies; i ++) {
            int[] values = getAllDataForSpecies(tree.getNode(i).getID());
            System.arraycopy(values, 0, arr, i * nrOfTraits, nrOfTraits);
        }
    }

    // getters
    public int getSpeciesNr () { return nrOfSpecies; }

    public int getTraitNr () { return nrOfTraits; }

    // get traitIdx-th trait values for speciesIdx-th species
    public int getTraitDataForSpecies(int speciesIdx, int traitIdx) {
        return discreteDataArr[speciesIdx * nrOfTraits + traitIdx];
    }

    protected void initiateLiabilities() {}

    public double[] getLiabilities () {
        return liabilitiesInput.get().getDoubleValues();
    }

    public int[] getDiscreteDataArr() {
        return discreteDataArr;
    }

    public double getLogP() { return logP; }

    @Override
    public List<String> getArguments() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public List<String> getConditions() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        // TODO Auto-generated method stub

    }

}
