package contraband.prunelikelihood;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.DataType;
import beast.evolution.tree.Tree;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

// This class implements mappings from discrete traits to continuous random variables.
public abstract class ThresholdModel extends Distribution {
    final public Input<RealParameter> liabilitiesInput = new Input<>("liability", "Continuous random variables underneath observed traits.", Input.Validate.REQUIRED);
    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Input.Validate.REQUIRED);

    private int[] discreteData;

    @Override
    public void initAndValidate() {

        Alignment data = dataInput.get();
        int nrOfStates = data.getDataType().getStateCount();
        int nrOfTraits = data.getSiteCount();
        List<Integer> n = data.getStateCounts();

        List<Sequence> seq = data.sequenceInput.get();
        Sequence seq1 = seq.get(0);validateInputs();
        DataType dy = data.getDataType();
        List<Integer> j = seq1.getSequence(dy);
        int nTaxa = data.getTaxonCount();
        Tree tree = treeInput.get();
        int nSpecies = tree.getLeafNodeCount();
        liabilitiesInput.get().setDimension(nrOfTraits * nSpecies);

        discreteData = new int[nrOfTraits * nSpecies];
        String[] speciesNames = new String[nSpecies];
        for (int i = 0; i < nSpecies; i ++) {
            speciesNames[i] = tree.getNode(i).getID();
        }

        System.out.println(Arrays.toString(speciesNames));


    }

    public double[] getLiabilities() {
        return liabilitiesInput.get().getDoubleValues();
    }



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
