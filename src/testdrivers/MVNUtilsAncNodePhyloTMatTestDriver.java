package testdrivers;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import beast.util.TreeParser;
import contraband.utils.GeneralUtils;
import contraband.math.MVNUtils;

public class MVNUtilsAncNodePhyloTMatTestDriver {

	public static void main(String[] args) {
		String treeStr = "((sp1:4.0,sp2:3.0):7.0,((sp3:2.0,sp4:1.0):5.0,sp5:8.0):6.0)";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		int nLeaves = myTree.getLeafNodeCount();
		int nAncNodes = myTree.getInternalNodeCount();
		RealMatrix phyloAncNodeWTMat = MatrixUtils.createRealMatrix(nAncNodes-1, nLeaves);
		
		MVNUtils.populateAncNodePhyloTMatrix(myTree, phyloAncNodeWTMat);
		GeneralUtils.displayRealMatrix(phyloAncNodeWTMat);
	}
}
