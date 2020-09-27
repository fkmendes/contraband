package testdrivers;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.coalescent.CoalCorrection;
import contraband.utils.GeneralUtils;

/**
 * @author Fabio K. Mendes
 */

public class CoalCorrectionTestDriver {

	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,(sp4:1.0,sp5:1.0):2.0);";
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);
		String[] spNamesInPhyloTMatOrder = new String[myTree.getLeafNodeCount()];

		// pop sizes
		Double[] popSizesInput = new Double[] { 0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.75, 1.0, 0.5 };
		RealParameter popSizes = new RealParameter(popSizesInput);
		
		CoalCorrection coal = new CoalCorrection();
		coal.initByName("tree", myTree, "popSizes", popSizes);
		double[][] correctedPhyloTMat = coal.getCorrectedPhyloTMat(spNamesInPhyloTMatOrder);
		
		GeneralUtils.display2DArray(correctedPhyloTMat);
	}
}
