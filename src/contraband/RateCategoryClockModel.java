package contraband;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * @author Alexei Drummond
 */
@Description("Local Strict Clock Model. There are nCat clock rate categories, and each branch is assigned to one of the nCat rate categories.")
public class RateCategoryClockModel extends BranchRateModel.Base {

    final public Input<IntegerParameter> rateCatAssignInput = new Input<>("rateCatAssign", "the rate category for each node in the tree.", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateValuesInput = new Input<>("rates", "the rate parameters associated with each category.", Input.Validate.REQUIRED);
    final public Input<Integer> nCatInput = new Input<>("nCat", "the number of rate categories.", Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);

    Tree tree;
    
    @Override
    public void initAndValidate() {
        tree = treeInput.get();

        IntegerParameter rateCategories = rateCatAssignInput.get();

        if (rateCategories.getDimension() != tree.getNodeCount()) {
            Log.warning.println("LocalStrictClockModel::Setting dimension of indicators to " + (tree.getNodeCount()));
            rateCategories.setDimension(tree.getNodeCount());
        }

        int nCat = nCatInput.get();
        for (int i = 0; i < rateCategories.getDimension(); i++) {
            if (rateCategories.getValue(i) > nCat) {
                Log.warning.println("Number of rate categories in assignment input exceeds the user-defined number of categories =" + nCat + " . Will fix that number.");
                nCat = rateCategories.getValue(i) + 1;
            }
        }

        RealParameter rates = rateValuesInput.get();

        /*
         * Note that this makes it necessary to set lower and upper to -Infinity and Infinity if
         * you expect that whatever "rate" represents can be a negative value.
         *
         * Otherwise, the lower will default to 0.0, and the logP of the prior of that parameter will
         * be -Infinity...
         */
        if (rates.lowerValueInput.get() == null || rates.lowerValueInput.get() < 0.0) {
            rates.setLower(0.0);
        }
        if (rates.upperValueInput.get() == null || rates.upperValueInput.get() < 0.0) {
            rates.setUpper(Double.MAX_VALUE);
        }
        if (rates.getDimension() != nCat) {
        	Log.warning.println("RandomLocalClockModel::Setting dimension of rates to " + nCat);
            rates.setDimension(nCat);
        }
    }

	/*
	 * Getters
	 */
    public int getNCat() {
    	return nCatInput.get();
    }
    
	public Double[] getColorValues() {
		return rateValuesInput.get().getValues();
	}
	
	public Integer[] getColorAssignments() {
		return rateCatAssignInput.get().getValues();
	}
	
    @Override
    public double getRateForBranch(Node node) {
        int index = rateCatAssignInput.get().getValue(node.getNr());
        return rateValuesInput.get().getArrayValue(index);
    }
}
