package testdrivers;

import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import beast.evolution.tree.Node;
import beast.util.TreeParser;
import contraband.GeneralUtils;
import contraband.MVNUtils;

public class RateMatrixTestDriver {

	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1[&Rate=1]:1.0, sp2[&Rate=1]:1.0)[&Rate=1]:1.0, sp3[&Rate=2.0]:2.0)[&Rate=2.0]:1.0, sp4[&Rate=3.0]:3.0)[&Rate=2.0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);

		int nodeCount = myTree.getNodeCount();
		int n = myTree.getLeafNodeCount();
		Double[] nodeRates = new Double[nodeCount];
		Node rootNode = myTree.getRoot();
		myTree.getMetaData(rootNode, nodeRates, "Rate");
		RealMatrix rateMatrix = new Array2DRowRealMatrix(n, n);
		
		MVNUtils.populateRateMatrix(myTree, nodeRates, rateMatrix);
		GeneralUtils.displayRealMatrix(rateMatrix);
//		for (int i=0; i<rateMatrix.length; ++i) {
//			System.out.println(Arrays.toString(rateMatrix[i]));
//		}
	}
}
