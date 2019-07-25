package contraband;


import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * @author Alexei Drummond
 */
@Description("Local Strict Clock Model. There are k clock rates, and each branch is assigned to one of the k rates.")
public class LocalStrictClockModel extends BranchRateModel.Base {

    final public Input<IntegerParameter> rateCategoryParamInput =
            new Input<>("rateCategories",
                    "the rate category for each node in the tree.",
                    Input.Validate.REQUIRED);
    final public Input<RealParameter> rateParamInput =
            new Input<>("rates",
                    "the rate parameters associated with each category.",
                    Input.Validate.REQUIRED);

    final public Input<Integer> kInput =
            new Input<>("k",
                    "the number of rate categories.",
                    Input.Validate.REQUIRED);

    final public Input<Tree> treeInput =
            new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);

    Tree tree;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();

        IntegerParameter rateCategories = rateCategoryParamInput.get();

        if (rateCategories.getDimension() != tree.getNodeCount()) {
            Log.warning.println("LocalStrictClockModel::Setting dimension of indicators to " + (tree.getNodeCount()));
            rateCategories.setDimension(tree.getNodeCount());
        }

        int k = kInput.get();
        for (int i = 0; i < rateCategories.getDimension(); i++) {
            if (rateCategories.getValue(i) >= k) {
                Log.warning.println("rate categories exceed the given k =" + k);
                k = rateCategories.getValue(i) + 1;

            }
        }

        RealParameter rates = rateParamInput.get();
        if (rates.lowerValueInput.get() == null || rates.lowerValueInput.get() < 0.0) {
            rates.setLower(0.0);
        }
        if (rates.upperValueInput.get() == null || rates.upperValueInput.get() < 0.0) {
            rates.setUpper(Double.MAX_VALUE);
        }
        if (rates.getDimension() != k) {
        	Log.warning.println("RandomLocalClockModel::Setting dimension of rates to " + k);
            rates.setDimension(k);
        }
    }

    @Override
    public double getRateForBranch(Node node) {
        int index = rateCategoryParamInput.get().getValue(node.getNr());
        return rateParamInput.get().getArrayValue(index);
    }
}
