package contraband;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * @author Alexei Drummond
 */
@Description("Random Local Color Model.")
@Citation(value =
        "Drummond AJ, Suchard MA (2010) Bayesian random local clocks, or one rate to rule them all. BMC biology 8, 114.",
        DOI = "10.1186/1741-7007-8-114",
        year = 2010, firstAuthorSurname = "drummond")
public class RandomLocalColorModel extends BranchRateModel.Base {

    final public Input<BooleanParameter> indicatorParamInput =
            new Input<>("indicators",
                    "the indicators associated with nodes in the tree for sampling of individual color changes among branches.",
                    Input.Validate.REQUIRED);
    final public Input<RealParameter> colorParamInput =
            new Input<>("colors",
                    "the color (rates, optima, etc.) parameters associated with nodes in the tree for sampling of individual colors among branches.",
                    Input.Validate.REQUIRED);

    final public Input<Tree> treeInput =
            new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);

    final public Input<Boolean> colorsAreMultipliersInput =
            new Input<>("colorsAreMultipliers", "true if the colors should be treated as multipliers (default false).",
                    false);

    final public Input<Boolean> scalingInput =
            new Input<>("scaling", "if false, then ignore meanColor input and leave colors unscaled.",
                    true, Input.Validate.OPTIONAL);

    final public Input<Boolean> includeRootInput =
            new Input<>("includeRoot", "if true, then the root can take on an arbitrary color, otherwise the root branch has color 1.0.",
                    false, Input.Validate.OPTIONAL);

    Tree tree;
    RealParameter meanColor;
    boolean scaling = true;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();

        scaling = scalingInput.get();

        BooleanParameter indicators = indicatorParamInput.get();

        int colorSize = tree.getNodeCount();
        int indicatorSize = tree.getNodeCount() - 1;

        if (!includeRootInput.get()) colorSize -= 1;

        if (indicators.getDimension() != indicatorSize) {
            Log.warning("RandomLocalColorModel::Setting dimension of indicators to " + indicatorSize);
            indicators.setDimension(indicatorSize);
        }

        unscaledBranchColors = new double[tree.getNodeCount()];

        RealParameter colors = colorParamInput.get();

        // Note: no restrictions for color values; rates in the
        // original relaxed clock model can only be > 0.0, but
        // here we allow negative colors, as often people do things
        // in log scale

        if (colors.getDimension() != colorSize) {
        	Log.warning("RandomLocalClockModel::Setting dimension of colors to " + colorSize);
            colors.setDimension(colorSize);
        }

        colorsAreMultipliers = colorsAreMultipliersInput.get();

        meanColor = meanColorInput.get();
        if (meanColor == null) {
            meanColor = new RealParameter("1.0");
        }
    }

    /**
     * This is a recursive function that does the work of
     * calculating the unscaled branch colors across the tree
     * taking into account the indicator variables.
     *
     * @param node the node
     * @param color the color of the parent node
     */
    private void calculateUnscaledBranchColors(Node node, double color, BooleanParameter indicators, RealParameter colors) {

        int nodeNumber = getNr(node);

        if (!node.isRoot()) {
            if (indicators.getValue(nodeNumber)) {
                if (colorsAreMultipliers) {
                    color *= colors.getValue(nodeNumber);
                } else {
                    color = colors.getValue(nodeNumber);
                }
            }
        }
        unscaledBranchColors[nodeNumber] = color;

        if (!node.isLeaf()) {
            calculateUnscaledBranchColors(node.getLeft(), color, indicators, colors);
            calculateUnscaledBranchColors(node.getRight(), color, indicators, colors);
        }
    }

    private void recalculateScaleFactor() {
        BooleanParameter indicators = indicatorParamInput.get();
        RealParameter colors = colorParamInput.get();

        double rootColor = 1.0;
        if (includeRootInput.get()) rootColor = colors.getValue(tree.getRoot().getNr());

        calculateUnscaledBranchColors(tree.getRoot(), rootColor, indicators, colors);

        if (scaling) {

            double timeTotal = 0.0;
            double branchTotal = 0.0;

            for (int i = 0; i < tree.getNodeCount(); i++) {
                Node node = tree.getNode(i);
                if (!node.isRoot()) {

                    double branchInTime = node.getParent().getHeight() - node.getHeight();

                    double branchLength = branchInTime * unscaledBranchColors[node.getNr()];

                    timeTotal += branchInTime;
                    branchTotal += branchLength;
                }
            }

            scaleFactor = timeTotal / branchTotal;

            scaleFactor *= meanColor.getValue();
        } else {
            scaleFactor = 1.0;
        }
    }

    @Override
    public double getRateForBranch(Node node) {
        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads
    	synchronized (this) {
    		if (recompute) {
                recalculateScaleFactor();
                recompute = false;
			}
        }

        return unscaledBranchColors[getNr(node)] * scaleFactor;
    }

    private int getNr(Node node) {
        int nodeNr = node.getNr();
        if (nodeNr > tree.getRoot().getNr()) {
            nodeNr--;
        }
        return nodeNr;
    }

    @Override
    protected boolean requiresRecalculation() {
        // this is only called if any of its inputs is dirty, hence we need to recompute
        recompute = true;
        return true;
    }

    @Override
    protected void store() {
        recompute = true;
        super.store();
    }

    @Override
    protected void restore() {
        recompute = true;
        super.restore();
    }

    private boolean recompute = true;
    double[] unscaledBranchColors;
    double scaleFactor = 1.0;
    boolean colorsAreMultipliers = false;
}
