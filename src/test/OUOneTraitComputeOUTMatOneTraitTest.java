package test;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import beast.util.TreeParser;

import java.util.ArrayList;
import java.util.List;

import contraband.GeneralUtils;

import beast.evolution.tree.Node;
import contraband.MVNUtils;
import contraband.OUUtils;

public class OUOneTraitComputeOUTMatOneTraitTest {
	
	final static double EPSILON = 1e-4;
	
//	For the covariance matrices we can assume two different hipothesis:
//		    'F' suffix refers to assuming a fixed parameter Root
//		    'R' suffix refers to assuming a random variable Root (stationary distribution)
//	For the weight matrices we can assume two different hipothesis:
//		    'I' suffix refers to isolating the root optimum value in the weight Matrix
//		    'M' suffix refers to merging the root parameter with the optimum parameter associated with the eldest selective regime
	
// 	Every test has four different outputs according to the combination of the previous situations: FI, FM, RI, RM	
//	 1, 2, 3 refers to Test 1, Test 2 and Test 3 respectively
	
	private static double[][] ouTMat1FI, ouTMat1FM, ouTMat1RI, ouTMat1RM;
	private static double[][] ouTMat2FI, ouTMat2FM, ouTMat2RI, ouTMat2RM;
	private static double[][] ouTMat3FI, ouTMat3FM, ouTMat3RI, ouTMat3RM;
	
	/*
	 * Test 1 contains 4 asserts. It checks the OU covariance matrix values in 4 different scenarios.
	 * 
	 * (1) Fixed root trait value, root optimum separate from the "primary optimum" (the primary optimum is
	 * the "ancestral" optimum that existed before the root, and kept being the optimum down one or both of 
	 * the branches coming from the root)
	 * 
	 * (2) Fixed root trait value, root optimum is set to be the same as whatever the primary optimum is.
	 * The primary optimum is passed by the user as the optimum of the root (some regimem is picked and that
	 * is both the primary optimum and the root optimum).
	 * 	  
	 * NOTE: Fixed root trait values means the ancestral state of the whole tree (the ancestral state at the root)
	 * is a parameter that we are sampling from during MCMC. That's the "normal" MVN likelihood as described by
	 * Felsenstein -- the x_0 parameter.
	 * 
	 * (3) Random root trait value, root optimum separate from the primary optimum
	 * 
	 * (4) Random root trait value, root optimum is set to be the same as whatever the primary optimum is (which is
	 * defined by the user through the root meta data in the tree string).
	 * 
	 * NOTE: Random root trait value means the ancestral state of the whole tree (the ancestral state at the root)
	 * is not a parameter that we sample from during MCMC. Instead, it comes from some stationary distribution.
	 * The density function of that distribution is worked into the covariance matrix of OU (and ultimately into
	 * the likelihood function of OU). That is why the covariance formulas for (F)ixed and (R)andom root cases
	 * are different. 
	 */
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
		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
		
	    // OU T matrix that will store the results (which will be asserted)
		ouTMat1FI = new double[n][n];
		ouTMat1FM = new double[n][n];
		ouTMat1RI = new double[n][n];
		ouTMat1RM = new double[n][n];
	
		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMat1FI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMat1FM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMat1RI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMat1RM, false);
		
		// Let's print stuff
		GeneralUtils.scalarBy2DArray(ouTMat1FI, sigsqFI);
		GeneralUtils.scalarBy2DArray(ouTMat1FM, sigsqFM);
		GeneralUtils.scalarBy2DArray(ouTMat1RI, sigsqRI);
		GeneralUtils.scalarBy2DArray(ouTMat1RM, sigsqRM);
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
		
		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
		
		// OU T matrix that will store the results (which will be asserted)
		ouTMat2FI = new double[n][n];
		ouTMat2FM = new double[n][n];
		ouTMat2RI = new double[n][n];
		ouTMat2RM = new double[n][n];
		
		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMat2FI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMat2FM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMat2RI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMat2RM, false);
		
		// Let's print stuff
		GeneralUtils.scalarBy2DArray(ouTMat2FI, sigsqFI);
		GeneralUtils.scalarBy2DArray(ouTMat2FM, sigsqFM);
		GeneralUtils.scalarBy2DArray(ouTMat2RI, sigsqRI);
		GeneralUtils.scalarBy2DArray(ouTMat2RM, sigsqRM);
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
		
		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
		
		// OU T matrix that will store the results (which will be asserted)
		ouTMat3FI = new double[n][n];
		ouTMat3FM = new double[n][n];
		ouTMat3RI = new double[n][n];
		ouTMat3RM = new double[n][n];
		
		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMat3FI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMat3FM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMat3RI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMat3RM, false);
		
		// Let's print stuff
		GeneralUtils.scalarBy2DArray(ouTMat3FI, sigsqFI);
		GeneralUtils.scalarBy2DArray(ouTMat3FM, sigsqFM);
		GeneralUtils.scalarBy2DArray(ouTMat3RI, sigsqRI);
		GeneralUtils.scalarBy2DArray(ouTMat3RM, sigsqRM);
	}

	@Test	// Test 1
	public void againstRvarOUTest1 () {
		
		Assert.assertEquals(1.563088e-29, 	ouTMat1FI[0][1], EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouTMat1FI[2][2], EPSILON);
		Assert.assertEquals(1.221620e-56, 	ouTMat1FI[1][2], EPSILON);
		Assert.assertEquals(1.563088e-29, 	ouTMat1FI[1][0], EPSILON);
		
		Assert.assertEquals(1.563088e-29, 	ouTMat1FM[0][1], EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouTMat1FM[2][2], EPSILON);
		Assert.assertEquals(1.221620e-56, 	ouTMat1FM[1][2], EPSILON);
		Assert.assertEquals(1.563088e-29, 	ouTMat1FM[1][0], EPSILON);
		
		Assert.assertEquals(4.417827e-40, 	ouTMat1RI[0][1], EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouTMat1RI[2][2], EPSILON);
		Assert.assertEquals(9.758589e-78, 	ouTMat1RI[1][2], EPSILON);
		Assert.assertEquals(4.417827e-40, 	ouTMat1RI[1][0], EPSILON);
		
		Assert.assertEquals(4.417827e-40, 	ouTMat1RM[0][1], EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouTMat1RM[2][2], EPSILON);
		Assert.assertEquals(9.758589e-78, 	ouTMat1RM[1][2], EPSILON);
		Assert.assertEquals(4.417827e-40, 	ouTMat1RM[1][0], EPSILON);
	}
	
	@Test	// Test 2
	public void againstRvarOUTest2 () {
		
		Assert.assertEquals(0.2711105, 		ouTMat2FI[0][1], EPSILON);
		Assert.assertEquals(2.67973939, 	ouTMat2FI[2][2], EPSILON);
		Assert.assertEquals(0.02357370, 	ouTMat2FI[1][2], EPSILON);
		Assert.assertEquals(0.27111053, 	ouTMat2FI[1][0], EPSILON);
		
		Assert.assertEquals(2.475727, 		ouTMat2FM[0][1], EPSILON);
		Assert.assertEquals(6.189318, 		ouTMat2FM[2][2], EPSILON);
		Assert.assertEquals(1.237863, 		ouTMat2FM[1][2], EPSILON);
		Assert.assertEquals(2.475727, 		ouTMat2FM[1][0], EPSILON);
		
		Assert.assertEquals(0.20187480, 	ouTMat2RI[0][1], EPSILON);
		Assert.assertEquals(2.672040045, 	ouTMat2RI[2][2], EPSILON);
		Assert.assertEquals(0.015251805, 	ouTMat2RI[1][2], EPSILON);
		Assert.assertEquals(0.201874800, 	ouTMat2RI[1][0], EPSILON);
		
		Assert.assertEquals(0.57628838, 	ouTMat2RM[0][1], EPSILON);
		Assert.assertEquals(4.82816083, 	ouTMat2RM[2][2], EPSILON);
		Assert.assertEquals(0.06878567, 	ouTMat2RM[1][2], EPSILON);
		Assert.assertEquals(0.57628838, 	ouTMat2RM[1][0], EPSILON);
	}
	
	@Test 	// Test 3
	public void againstRvarOUTest3 () {
		
		Assert.assertEquals(4.640928e-07, 	ouTMat3FI[0][1], EPSILON);
		Assert.assertEquals(7.428888e-01, 	ouTMat3FI[2][2], EPSILON);
		Assert.assertEquals(3.668135e-10, 	ouTMat3FI[1][2], EPSILON);
		Assert.assertEquals(4.640928e-07, 	ouTMat3FI[1][0], EPSILON);
		
		Assert.assertEquals(4.640928e-07, 	ouTMat3FM[0][1], EPSILON);
		Assert.assertEquals(7.428888e-01, 	ouTMat3FM[2][2], EPSILON);
		Assert.assertEquals(3.668135e-10, 	ouTMat3FM[1][2], EPSILON);
		Assert.assertEquals(4.640928e-07, 	ouTMat3FM[1][0], EPSILON);
		
		Assert.assertEquals(1.139565e-07, 	ouTMat3RI[0][1], EPSILON);
		Assert.assertEquals(7.428742e-01, 	ouTMat3RI[2][2], EPSILON);
		Assert.assertEquals(4.463248e-11, 	ouTMat3RI[1][2], EPSILON);
		Assert.assertEquals(1.139565e-07, 	ouTMat3RI[1][0], EPSILON);
		
		Assert.assertEquals(1.139565e-07, 	ouTMat3RM[0][1], EPSILON);
		Assert.assertEquals(7.428742e-01, 	ouTMat3RM[2][2], EPSILON);
		Assert.assertEquals(4.463248e-11, 	ouTMat3RM[1][2], EPSILON);
		Assert.assertEquals(1.139565e-07, 	ouTMat3RM[1][0], EPSILON);
	}
}
