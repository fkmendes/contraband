package contraband.test;

import beast.base.evolution.tree.Tree;
import contraband.math.MatrixUtilsContra;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;

import beast.base.evolution.tree.TreeParser;

import java.util.ArrayList;
import java.util.List;

import beast.base.evolution.tree.Node;
import contraband.utils.MVNUtils;
import contraband.utils.OUUtils;

/**
 * @author Pau Bravo
 */

/*
 * This class contains unit tests for OUUtils
 */
public class OUUtilsComputeOUTMatOneTraitTest {
	
	final static double EPSILON = 1e-4;

	/* 
	 * The matrix we are computing here (after multiplying by sigma^2)
	 * is what is referred to in Butler & King as the V (variance) matrix of
	 * the OU process.
	 * 
	 * We call it the T matrix here because initially we grab the evolutionary
	 * distances (without sigma^2), and then multiplied by sigma^2 in place --
	 * at this point, we are looking at the OU process V matrix.
	 * 
	 * R: rootIsRandVar=true
	 * F: rootIsRandVar=false (assumes equilibrium distribution at root)
	 * I: useRootMetaData=true
	 * M: useRootMetaData=false (one fewer parameter)
	 */
	private static RealMatrix ouTMat1FI, ouTMat1FM, ouTMat1RI, ouTMat1RM;
	private static RealMatrix ouTMat2FI, ouTMat2FM, ouTMat2RI, ouTMat2RM;
	private static RealMatrix ouTMat3FI, ouTMat3FM, ouTMat3RI, ouTMat3RM;

	/*
	 * (1) Small ultrametric tree, in all combinations of R, F, I and M
	 */
	@Test
	public void testComputeOUTMatSmallTree() {

		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// The values below are MLEs from mvMORPH
		double alphaFI = 31.20814;
		double alphaFM = 31.20814;
		double alphaRI = 43.35287;
		double alphaRM = 43.35287;
		
		double sigsqFI = 1.248328;		
		double sigsqFM = 1.248328;
		double sigsqRI = 1.734117;
		double sigsqRM = 1.734117;
		
		// Creating T matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		String[] spOrderInTMat = new String[n];
		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves, spOrderInTMat);
		
	    // OU T matrix that will store the results (which will be asserted)
		ouTMat1FI = new Array2DRowRealMatrix(n, n);
		ouTMat1FM = new Array2DRowRealMatrix(n, n);
		ouTMat1RI = new Array2DRowRealMatrix(n, n);
		ouTMat1RM = new Array2DRowRealMatrix(n, n);
	
		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMat1FI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMat1FM, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMat1RI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMat1RM, true);

		// multiplying by sigma^2 in place
		MatrixUtilsContra.scalarByRealMatrix(ouTMat1FI, sigsqFI);
		MatrixUtilsContra.scalarByRealMatrix(ouTMat1FM, sigsqFM);
		MatrixUtilsContra.scalarByRealMatrix(ouTMat1RI, sigsqRI);
		MatrixUtilsContra.scalarByRealMatrix(ouTMat1RM, sigsqRM);

		Assert.assertEquals(1.563088e-29, ouTMat1FI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(2.000003e-02, ouTMat1FI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(1.221620e-56, ouTMat1FI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.563088e-29, ouTMat1FI.getEntry(1, 0), EPSILON);

		Assert.assertEquals(1.563088e-29, ouTMat1FM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(2.000003e-02, ouTMat1FM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(1.221620e-56, ouTMat1FM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.563088e-29, ouTMat1FM.getEntry(1, 0), EPSILON);

		Assert.assertEquals(4.417827e-40, ouTMat1RI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(2.000003e-02, ouTMat1RI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(9.758589e-78, ouTMat1RI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(4.417827e-40, ouTMat1RI.getEntry(1, 0), EPSILON);
	}

	/*
	 * (2) Small non-ultrametric tree, in all combinations of R, F, I and M
	 */
	@Test
	public void testComputeOUTMatSmallTreeNonUltra() {

		String treeStrNonUltra = "(((sp1:2.0,sp2:1.0):1.0,sp3:4.0):1.0,sp4:3.0);";
		Tree myTree = new TreeParser(treeStrNonUltra, false, false, true, 0);
		
		// The values below are MLEs from mvMORPH
		double alphaFI = 0.7465763;
		double alphaFM = 1.40338e-08;
		double alphaRI = 0.8609833;
		double alphaRM = 0.7085376;
		
		double sigsqFI = 4.003551;				
		double sigsqFM = 1.237864;
		double sigsqRI = 4.601164;
		double sigsqRM = 6.841867;
		
		// Creating T matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];	
		double[] nodeToRootPaths = new double[myTree.getNodeCount()]; 
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		String[] spOrderInTMat = new String[n];
		
		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves, spOrderInTMat);
		
		// OU T matrix that will store the results (which will be asserted)
		ouTMat2FI = new Array2DRowRealMatrix(n, n);
		ouTMat2FM = new Array2DRowRealMatrix(n, n);
		ouTMat2RI = new Array2DRowRealMatrix(n, n);
		ouTMat2RM = new Array2DRowRealMatrix(n, n);
		
		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMat2FI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMat2FM, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMat2RI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMat2RM, true);

		// multiplying by sigma^2 in place
		MatrixUtilsContra.scalarByRealMatrix(ouTMat2FI, sigsqFI);
		MatrixUtilsContra.scalarByRealMatrix(ouTMat2FM, sigsqFM);
		MatrixUtilsContra.scalarByRealMatrix(ouTMat2RI, sigsqRI);
		MatrixUtilsContra.scalarByRealMatrix(ouTMat2RM, sigsqRM);

		Assert.assertEquals(0.2711105, ouTMat2FI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(2.67973939, ouTMat2FI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0.02357370, ouTMat2FI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(0.27111053, ouTMat2FI.getEntry(1, 0), EPSILON);

		Assert.assertEquals(2.475727, ouTMat2FM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(6.189318, ouTMat2FM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(1.237863, ouTMat2FM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(2.475727, ouTMat2FM.getEntry(1, 0), EPSILON);

		Assert.assertEquals(0.20187480, ouTMat2RI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(2.672040045, ouTMat2RI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0.015251805, ouTMat2RI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(0.201874800, ouTMat2RI.getEntry(1, 0), EPSILON);

		Assert.assertEquals(0.57628838, ouTMat2RM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(4.82816083, ouTMat2RM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0.06878567, ouTMat2RM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(0.57628838, ouTMat2RM.getEntry(1, 0), EPSILON);
	}
	
	/*
	 * (3) Large non-ultrametric tree, in all combinations of R, F, I and M
	 */
	@Test
	public void testComputeOUTMatLargeTreeNonUltra() {

		String treeStrNonUltraBig = "(((((sp1:1.0,sp2:1.0):1.0,sp3:1.0):2.0,(sp4:1.0,sp5:1.0):3.0):2.0, sp6:6.0):1.0,sp7:3.0);";
		Tree myTree = new TreeParser(treeStrNonUltraBig, false, false, true, 0);

		// The values below are MLEs from mvMORPH
		double alphaFI = 7.142986;
		double alphaFM = 7.142986;
		double alphaRI = 7.84511;
		double alphaRM = 7.84511;

		double sigsqFI = 10.61289;
		double sigsqFM = 10.61289;
		double sigsqRI = 11.65586;
		double sigsqRM = 11.65586;

		// Creating T matrix of the tree
		int n = myTree.getLeafNodeCount(); // Delete later
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		String[] spOrderInTMat = new String[n];

		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves, spOrderInTMat);

		// OU T matrix that will store the results (which will be asserted)
		ouTMat3FI = new Array2DRowRealMatrix(n, n);
		ouTMat3FM = new Array2DRowRealMatrix(n, n);
		ouTMat3RI = new Array2DRowRealMatrix(n, n);
		ouTMat3RM = new Array2DRowRealMatrix(n, n);

		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMat3FI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMat3FM, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMat3RI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMat3RM, true);

		// multiplying by sigma^2 in place
		MatrixUtilsContra.scalarByRealMatrix(ouTMat3FI, sigsqFI);
		MatrixUtilsContra.scalarByRealMatrix(ouTMat3FM, sigsqFM);
		MatrixUtilsContra.scalarByRealMatrix(ouTMat3RI, sigsqRI);
		MatrixUtilsContra.scalarByRealMatrix(ouTMat3RM, sigsqRM);

		Assert.assertEquals(4.640928e-07, ouTMat3FI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(7.428888e-01, ouTMat3FI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(3.668135e-10, ouTMat3FI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(4.640928e-07, ouTMat3FI.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(4.640928e-07, ouTMat3FM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(7.428888e-01, ouTMat3FM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(3.668135e-10, ouTMat3FM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(4.640928e-07, ouTMat3FM.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(1.139565e-07, ouTMat3RI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(7.428742e-01, ouTMat3RI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(4.463248e-11, ouTMat3RI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.139565e-07, ouTMat3RI.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(1.139565e-07, ouTMat3RM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(7.428742e-01, ouTMat3RM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(4.463248e-11, ouTMat3RM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.139565e-07, ouTMat3RM.getEntry(1, 0), EPSILON);
	}
}
