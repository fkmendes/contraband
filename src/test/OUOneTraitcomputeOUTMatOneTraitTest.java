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

public class OUOneTraitcomputeOUTMatOneTraitTest {
	
	final static double EPSILON = 1e-4;
	
//	For the covariance matrices we can assume two different hipothesis:
//		    'F' suffix refers to assuming a fixed parameter Root
//		    'R' suffix refers to assuming a random variable Root (stationary distribution)
//	For the weight matrices we can assume two different hipothesis:
//		    'I' suffix refers to isolating the root optimum value in the weight Matrix
//		    'M' suffix refers to merging the root parameter with the optimum parameter associated with the eldest selective regime
	
// 	Every test has four different outputs according to the combination of the previous situations: FI, FM, RI, RM	
//	 1, 2, 3 refers to Test 1, Test 2 and Test 3 respectively
	
	private static double[][] ouCov1FI, ouCov1FM, ouCov1RI, ouCov1RM;
	private static double[][] ouCov2FI, ouCov2FM, ouCov2RI, ouCov2RM;
	private static double[][] ouCov3FI, ouCov3FM, ouCov3RI, ouCov3RM;
	
	
	@BeforeClass	// Test1: Ultrametric tree (3 regimes)
	public static void setUPTest1() throws Exception {
		
	
		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		double alphaFI = 31.20814;				// mvMORPH values
		double alphaFM = 31.20814;
		double alphaRI = 43.35287;
		double alphaRM = 43.35287;
		
		double sigsqFI = 1.248328;		
		double sigsqFM = 1.248328;
		double sigsqRI = 1.734117;
		double sigsqRM = 1.734117;
		
			// Creating phylogenetic covariance matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		
		MVNUtils.populateVcvMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
		
			// OU covariance matrices
		ouCov1FI = new double[n][n];
		ouCov1FM = new double[n][n];
		ouCov1RI = new double[n][n];
		ouCov1RM = new double[n][n];
	
		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouCov1FI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouCov1FM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouCov1RI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouCov1RM, false);
		
		GeneralUtils.scalarBy2DArray(ouCov1FI, sigsqFI);
		GeneralUtils.scalarBy2DArray(ouCov1FM, sigsqFM);
		GeneralUtils.scalarBy2DArray(ouCov1RI, sigsqRI);
		GeneralUtils.scalarBy2DArray(ouCov1RM, sigsqRM);
		
	}
	
	@BeforeClass	// Test2: non-ultrametric tree (1 regimes)
	public static void setUPTest2() throws Exception {
		
		
		String treeStrNonUltra = "(((sp1[&Regime=0]:2.0, sp2[&Regime=0]:1.0)[&Regime=0]:1.0, sp3[&Regime=0]:4.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStrNonUltra, false, false, true, 0);
		
		double alphaFI = 0.7465763 ;				// mvMORPH values
		double alphaFM = 1.40338e-08;
		double alphaRI = 0.8609833;
		double alphaRM = 0.7085376;
		
		double sigsqFI = 4.003551;				
		double sigsqFM = 1.237864;
		double sigsqRI = 4.601164;
		double sigsqRM = 6.841867;
		
			// Creating phylogenetic covariance matrix of the tree
		int n = myTree.getLeafNodeCount(); // Delete later
		double[][] tMatInput = new double[n][n];	
		double[] nodeToRootPaths = new double[myTree.getNodeCount()]; 
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		
		MVNUtils.populateVcvMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
		
		// OU covariance matrices
		ouCov2FI = new double[n][n];
		ouCov2FM = new double[n][n];
		ouCov2RI = new double[n][n];
		ouCov2RM = new double[n][n];
		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouCov2FI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouCov2FM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouCov2RI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouCov2RM, false);
		
		GeneralUtils.scalarBy2DArray(ouCov2FI, sigsqFI);
		GeneralUtils.scalarBy2DArray(ouCov2FM, sigsqFM);
		GeneralUtils.scalarBy2DArray(ouCov2RI, sigsqRI);
		GeneralUtils.scalarBy2DArray(ouCov2RM, sigsqRM);
	}
	
	@BeforeClass	// Test3: non-ultrametric tree (5 regimes)
	public static void setUPTest3() throws Exception {
		
		String treeStrNonUltraBig = "(((((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=0]:1.0)[&Regime=0]:2.0, (sp4[&Regime=2]:1.0, sp5[&Regime=2]:1.0)[&Regime=2]:3.0)[&Regime=0]:2.0, sp6[&Regime=3]:6.0)[&Regime=0]:1.0, sp7[&Regime=4]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStrNonUltraBig, false, false, true, 0);
		
		double alphaFI = 7.142986;				// mvMORPH values
		double alphaFM = 7.142986;
		double alphaRI = 7.84511;
		double alphaRM = 7.84511;
		
		double sigsqFI = 10.61289;	
		double sigsqFM = 10.61289;
		double sigsqRI = 11.65586;
		double sigsqRM = 11.65586;
		
			// Creating phylogenetic covariance matrix of the tree
		int n = myTree.getLeafNodeCount(); // Delete later
		double[][] tMatInput = new double[n][n];	
		double[] nodeToRootPaths = new double[myTree.getNodeCount()]; 
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		
		MVNUtils.populateVcvMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
		
			// OU covariance matrices
		// OU covariance matrices
		ouCov3FI = new double[n][n];
		ouCov3FM = new double[n][n];
		ouCov3RI = new double[n][n];
		ouCov3RM = new double[n][n];
		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouCov3FI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouCov3FM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouCov3RI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouCov3RM, false);
		
		GeneralUtils.scalarBy2DArray(ouCov3FI, sigsqFI);
		GeneralUtils.scalarBy2DArray(ouCov3FM, sigsqFM);
		GeneralUtils.scalarBy2DArray(ouCov3RI, sigsqRI);
		GeneralUtils.scalarBy2DArray(ouCov3RM, sigsqRM);
	}

	@Test	// Test 1
	public void againstRvarOUTest1 () {
		
		Assert.assertEquals(1.563088e-29, 	ouCov1FI[0][1], EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouCov1FI[2][2], EPSILON);
		Assert.assertEquals(1.221620e-56, 	ouCov1FI[1][2], EPSILON);
		Assert.assertEquals(1.563088e-29, 	ouCov1FI[1][0], EPSILON);
		
		Assert.assertEquals(1.563088e-29, 	ouCov1FM[0][1], EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouCov1FM[2][2], EPSILON);
		Assert.assertEquals(1.221620e-56, 	ouCov1FM[1][2], EPSILON);
		Assert.assertEquals(1.563088e-29, 	ouCov1FM[1][0], EPSILON);
		
		Assert.assertEquals(4.417827e-40, 	ouCov1RI[0][1], EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouCov1RI[2][2], EPSILON);
		Assert.assertEquals(9.758589e-78, 	ouCov1RI[1][2], EPSILON);
		Assert.assertEquals(4.417827e-40, 	ouCov1RI[1][0], EPSILON);
		
		Assert.assertEquals(4.417827e-40, 	ouCov1RM[0][1], EPSILON);
		Assert.assertEquals(2.000003e-02, 	ouCov1RM[2][2], EPSILON);
		Assert.assertEquals(9.758589e-78, 	ouCov1RM[1][2], EPSILON);
		Assert.assertEquals(4.417827e-40, 	ouCov1RM[1][0], EPSILON);
	}
	
	@Test	// Test 2
	public void againstRvarOUTest2 () {
		
		Assert.assertEquals(0.2711105, 		ouCov2FI[0][1], EPSILON);
		Assert.assertEquals(2.67973939, 	ouCov2FI[2][2], EPSILON);
		Assert.assertEquals(0.02357370, 	ouCov2FI[1][2], EPSILON);
		Assert.assertEquals(0.27111053, 	ouCov2FI[1][0], EPSILON);
		
		Assert.assertEquals(2.475727, 		ouCov2FM[0][1], EPSILON);
		Assert.assertEquals(6.189318, 		ouCov2FM[2][2], EPSILON);
		Assert.assertEquals(1.237863, 		ouCov2FM[1][2], EPSILON);
		Assert.assertEquals(2.475727, 		ouCov2FM[1][0], EPSILON);
		
		Assert.assertEquals(0.20187480, 	ouCov2RI[0][1], EPSILON);
		Assert.assertEquals(2.672040045, 	ouCov2RI[2][2], EPSILON);
		Assert.assertEquals(0.015251805, 	ouCov2RI[1][2], EPSILON);
		Assert.assertEquals(0.201874800, 	ouCov2RI[1][0], EPSILON);
		
		Assert.assertEquals(0.57628838, 	ouCov2RM[0][1], EPSILON);
		Assert.assertEquals(4.82816083, 	ouCov2RM[2][2], EPSILON);
		Assert.assertEquals(0.06878567, 	ouCov2RM[1][2], EPSILON);
		Assert.assertEquals(0.57628838, 	ouCov2RM[1][0], EPSILON);
	}
	
	@Test 	// Test 3
	public void againstRvarOUTest3 () {
		
		Assert.assertEquals(4.640928e-07, 	ouCov3FI[0][1], EPSILON);
		Assert.assertEquals(7.428888e-01, 	ouCov3FI[2][2], EPSILON);
		Assert.assertEquals(3.668135e-10, 	ouCov3FI[1][2], EPSILON);
		Assert.assertEquals(4.640928e-07, 	ouCov3FI[1][0], EPSILON);
		
		Assert.assertEquals(4.640928e-07, 	ouCov3FM[0][1], EPSILON);
		Assert.assertEquals(7.428888e-01, 	ouCov3FM[2][2], EPSILON);
		Assert.assertEquals(3.668135e-10, 	ouCov3FM[1][2], EPSILON);
		Assert.assertEquals(4.640928e-07, 	ouCov3FM[1][0], EPSILON);
		
		Assert.assertEquals(1.139565e-07, 	ouCov3RI[0][1], EPSILON);
		Assert.assertEquals(7.428742e-01, 	ouCov3RI[2][2], EPSILON);
		Assert.assertEquals(4.463248e-11, 	ouCov3RI[1][2], EPSILON);
		Assert.assertEquals(1.139565e-07, 	ouCov3RI[1][0], EPSILON);
		
		Assert.assertEquals(1.139565e-07, 	ouCov3RM[0][1], EPSILON);
		Assert.assertEquals(7.428742e-01, 	ouCov3RM[2][2], EPSILON);
		Assert.assertEquals(4.463248e-11, 	ouCov3RM[1][2], EPSILON);
		Assert.assertEquals(1.139565e-07, 	ouCov3RM[1][0], EPSILON);
	}

}
