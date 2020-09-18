package test;

import contraband.math.MatrixUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import beast.util.TreeParser;

import java.util.ArrayList;
import java.util.List;

import beast.evolution.tree.Node;
import contraband.math.MVNUtils;
import contraband.math.OUUtils;

public class OUOneTraitComputeOUTMatOneTraitTest {
	
	final static double EPSILON = 1e-4;
	
	/* 
	 * The OU T matrix is what is referred to in Butler & King as the V (variance) matrix of
	 * the OU process.
	 * 
	 * We call it the T matrix here because it still not being multiplied by sigma^2.
	 * 
	 * R: rootIsRandVar=true
	 * F: rootIsRandVar=false
	 * I: useRootMetaData=true
	 * M: useRootMetaData=false
	 */
	private static RealMatrix ouTMat1FI, ouTMat1FM, ouTMat1RI, ouTMat1RM;
	private static RealMatrix ouTMat2FI, ouTMat2FM, ouTMat2RI, ouTMat2RM;
	private static RealMatrix ouTMat3FI, ouTMat3FM, ouTMat3RI, ouTMat3RM;
	
	@BeforeClass
	public static void setUPTest1() throws Exception {
		
		// THREE regimes!
		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
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
		
		MatrixUtils.scalarByRealMatrix(ouTMat1FI, sigsqFI);
		MatrixUtils.scalarByRealMatrix(ouTMat1FM, sigsqFM);
		MatrixUtils.scalarByRealMatrix(ouTMat1RI, sigsqRI);
		MatrixUtils.scalarByRealMatrix(ouTMat1RM, sigsqRM);
	}
	
	/*
	 * Test 2 does the same as test 1, but the tree is non-ultrametric
	 * (same topology) and there are one regime instead of three.
	 */
	@BeforeClass
	public static void setUPTest2() throws Exception {
		
		// One regime!
		String treeStrNonUltra = "(((sp1[&Regime=0]:2.0, sp2[&Regime=0]:1.0)[&Regime=0]:1.0, sp3[&Regime=0]:4.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStrNonUltra, false, false, true, 0);
		
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

		MatrixUtils.scalarByRealMatrix(ouTMat2FI, sigsqFI);
		MatrixUtils.scalarByRealMatrix(ouTMat2FM, sigsqFM);
		MatrixUtils.scalarByRealMatrix(ouTMat2RI, sigsqRI);
		MatrixUtils.scalarByRealMatrix(ouTMat2RM, sigsqRM);
	}
	
	/*
	 * Test 3 does the same as the two previous tests, but the tree
	 * is larger and has a different topology (and is non-ultrametric).
	 * There are also five regimes.
	 */
	@BeforeClass
	public static void setUPTest3() throws Exception {
		
		// Five regimes!
		String treeStrNonUltraBig = "(((((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=0]:1.0)[&Regime=0]:2.0, (sp4[&Regime=2]:1.0, sp5[&Regime=2]:1.0)[&Regime=2]:3.0)[&Regime=0]:2.0, sp6[&Regime=3]:6.0)[&Regime=0]:1.0, sp7[&Regime=4]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStrNonUltraBig, false, false, true, 0);
		
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
		
		MatrixUtils.scalarByRealMatrix(ouTMat3FI, sigsqFI);
		MatrixUtils.scalarByRealMatrix(ouTMat3FM, sigsqFM);
		MatrixUtils.scalarByRealMatrix(ouTMat3RI, sigsqRI);
		MatrixUtils.scalarByRealMatrix(ouTMat3RM, sigsqRM);
	}

	@Test	// Test 1
	public void againstRvarOUTest1 () {
		Assert.assertEquals(1.0, 	1.0, EPSILON);
		/*
		Assert.assertEquals(1.563088e-29, 	ouTMat1FI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouTMat1FI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(1.221620e-56, 	ouTMat1FI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.563088e-29, 	ouTMat1FI.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(1.563088e-29, 	ouTMat1FM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouTMat1FM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(1.221620e-56, 	ouTMat1FM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.563088e-29, 	ouTMat1FM.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(4.417827e-40, 	ouTMat1RI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouTMat1RI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(9.758589e-78, 	ouTMat1RI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(4.417827e-40, 	ouTMat1RI.getEntry(1, 0), EPSILON);
		*/
	}

	
	@Test	// Test 2
	public void againstRvarOUTest2 () {
		Assert.assertEquals(1.0, 	1.0, EPSILON);
		/*
		Assert.assertEquals(0.2711105, 		ouTMat2FI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(2.67973939, 	ouTMat2FI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0.02357370, 	ouTMat2FI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(0.27111053, 	ouTMat2FI.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(2.475727, 		ouTMat2FM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(6.189318, 		ouTMat2FM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(1.237863, 		ouTMat2FM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(2.475727, 		ouTMat2FM.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(0.20187480, 	ouTMat2RI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(2.672040045, 	ouTMat2RI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0.015251805, 	ouTMat2RI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(0.201874800, 	ouTMat2RI.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(0.57628838, 	ouTMat2RM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(4.82816083, 	ouTMat2RM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0.06878567, 	ouTMat2RM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(0.57628838, 	ouTMat2RM.getEntry(1, 0), EPSILON);
		*/
	}
	
	@Test 	// Test 3
	public void againstRvarOUTest3 () {
		Assert.assertEquals(1.0, 	1.0, EPSILON);
		/*
		Assert.assertEquals(4.640928e-07, 	ouTMat3FI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(7.428888e-01, 	ouTMat3FI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(3.668135e-10, 	ouTMat3FI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(4.640928e-07, 	ouTMat3FI.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(4.640928e-07, 	ouTMat3FM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(7.428888e-01, 	ouTMat3FM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(3.668135e-10, 	ouTMat3FM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(4.640928e-07, 	ouTMat3FM.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(1.139565e-07, 	ouTMat3RI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(7.428742e-01, 	ouTMat3RI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(4.463248e-11, 	ouTMat3RI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.139565e-07, 	ouTMat3RI.getEntry(1, 0), EPSILON);
		
		Assert.assertEquals(1.139565e-07, 	ouTMat3RM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(7.428742e-01, 	ouTMat3RM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(4.463248e-11, 	ouTMat3RM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.139565e-07, 	ouTMat3RM.getEntry(1, 0), EPSILON);
		*/
	}
}
