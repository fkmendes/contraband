package testdrivers;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.clock.TreeToVCVMat;

/**
 * @author Fabio K. Mendes
 */

/*
 * Matches testRegimeManagerValues
 */
public class RegimeManagerTestDriver {

	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):2.0,(sp4:2.5,sp5:2.5):1.5);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
				
		// initializing data
		RealParameter colorValues = new RealParameter(new Double[] { 0.2, 0.4, 0.6, 0.8, 1.0 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 3, 3, 0, 0, 0, 2, 1, 4, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 5, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

		TreeToVCVMat colors = new TreeToVCVMat();
		colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);
		
		double[][] colorValuesMat = colors.getSpColorValuesMatOneTrait();
		for (int i=0; i<colorValuesMat.length; ++i) {
			for (int j=0; j<colorValuesMat[i].length; ++j) {
				String val = String.format("%.1f", colorValuesMat[i][j]);
				System.out.printf(val + " ");
			}

			System.out.println();
		}
		
		Node sp1Node = myTree.getNode(0);
		System.out.println(colors.getNodeColorValue(sp1Node, 0));
		
		Node sp3Node = myTree.getNode(2);
		System.out.println(colors.getNodeColorValue(sp3Node, 0));
	}
}
