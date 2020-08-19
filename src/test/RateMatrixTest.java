package test;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.evolution.tree.Node;
import beast.util.TreeParser;
import contraband.utils.GeneralUtils;
import contraband.math.MVNUtils;

public class RateMatrixTest {

	RealMatrix rateMatrix;
	double[] row1Expected, row1;
	double[] row2Expected, row2;
	double[] row3Expected, row3;
	double[] row4Expected, row4;
	
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((sp1[&Rate=1]:1.0, sp2[&Rate=1]:1.0)[&Rate=1]:1.0, sp3[&Rate=2.0]:2.0)[&Rate=2.0]:1.0, sp4[&Rate=3.0]:3.0)[&Rate=4.0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		int nodeCount = myTree.getNodeCount();
		int n = myTree.getLeafNodeCount();
		Double[] nodeRates = new Double[nodeCount];
		Node rootNode = myTree.getRoot();
		myTree.getMetaData(rootNode, nodeRates, "Rate");
		RealMatrix rateMatrix = new Array2DRowRealMatrix(n, n);
		
		MVNUtils.populateRateMatrix(myTree, nodeRates, rateMatrix);
		
		row1Expected = new double[] { 1.0, 1.0, 2.0, 4.0 }; // all elements should have weighted averages for their colors, this is wrong!
		row1 = rateMatrix.getRow(0);
		row2Expected = new double[] { 1.0, 1.0, 2.0, 4.0 };
		row2 = rateMatrix.getRow(1);
		row3Expected = new double[] { 2.0, 2.0, 2.0, 4.0 };
		row3 = rateMatrix.getRow(2);
		row4Expected = new double[] { 4.0, 4.0, 4.0, 3.0 };
		row4 = rateMatrix.getRow(3);
		
		GeneralUtils.displayRealMatrix(rateMatrix);
	}

	@Test
	public void test() {
		Assert.assertArrayEquals(row1Expected, row1, 0.0);
		Assert.assertArrayEquals(row2Expected, row2, 0.0);
		Assert.assertArrayEquals(row3Expected, row3, 0.0);
		Assert.assertArrayEquals(row4Expected, row4, 0.0);
	}

}
