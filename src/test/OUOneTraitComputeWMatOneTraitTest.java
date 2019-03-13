package test;

import java.util.List;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import beast.evolution.tree.Node;
import beast.util.TreeParser;

import contraband.OUUtils;

public class OUOneTraitComputeWMatOneTraitTest {

	final static double EPSILON = 1e-4;
	
	/* 
	 * F: root theta (=root mean=root ancestral trait value=root regime) is a separate parameter,
	 * we "F"ix the root value (condition OU on this value)
	 * R: root theta is not a parameter, it is a "R"andom value with a stationary distribution 
	 * that we integrate over
	 * 
	 * I: we set the root theta (=root mean=root ancestral trait value=root regime) to be its own "I"ndependent
	 * parameter, and estimate the
	 * M: we set the root theta to be the same ("M"erge) as the regime of one of its children, making it not be a parameter
	 */
	
	// We are using the same trees as in OUOneTraitcomputeOUTMatOneTraitTest.java)
	
	private static double[][] ouWeight1FI, ouWeight1FM, ouWeight1RI, ouWeight1RM;
	private static double[][] ouWeight2FI, ouWeight2FM, ouWeight2RI, ouWeight2RM;
	private static double[][] ouWeight3FI, ouWeight3FM, ouWeight3RI, ouWeight3RM;
	
	// Test 1: Ultrametric tree (3 regimes)
	@BeforeClass
	public static void setUPTest1() throws Exception {

		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();
		
		// The values below are MLEs from mvMORPH
		double alphaFI = 31.20814;
		double alphaFM = 31.20814;
		double alphaRI = 43.35287;
		double alphaRM = 43.35287;
		
		// Setting weight matrices dimensions
		int n = myTree.getLeafNodeCount();
		ouWeight1FI = new double[n][4];
		ouWeight1FM = new double[n][3];
		ouWeight1RI = new double[n][4];
		ouWeight1RM = new double[n][3];

	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 3, alphaFI, ouWeight1FI, false);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 3, alphaFM, ouWeight1FM, true);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 3, alphaRI, ouWeight1RI, false);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 3, alphaRM, ouWeight1RM, true);
	}
	
	// Test 2: Non-ultrametric tree (1 regimes)
	@BeforeClass
	public static void setUPTest2() throws Exception {
		
		String treeStrNonUltra = "(((sp1[&Regime=0]:2.0, sp2[&Regime=0]:1.0)[&Regime=0]:1.0, sp3[&Regime=0]:4.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStrNonUltra, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();

		// The values below are MLEs from mvMORPH
		double alphaFI = 0.7465763;
		double alphaFM = 1.40338e-08;
		double alphaRI = 0.8609833;
		double alphaRM = 0.7085376;
		
		// Setting weight matrices dimensions
		int n = myTree.getLeafNodeCount();
		ouWeight2FI = new double[n][2];
		ouWeight2FM = new double[n][1];
		ouWeight2RI = new double[n][2];
		ouWeight2RM = new double[n][1];

		OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 1, alphaFI, ouWeight2FI, false);
		OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 1, alphaFM, ouWeight2FM, true);
		OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 1, alphaRI, ouWeight2RI, false);
		OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 1, alphaRM, ouWeight2RM, true);
	}
	
	// Test 3: Non-ultrametric tree (5 regimes)
	@BeforeClass
	public static void setUPTest3() throws Exception {
		
		String treeStrNonUltraBig = "(((((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=0]:1.0)[&Regime=0]:2.0, (sp4[&Regime=2]:1.0, sp5[&Regime=2]:1.0)[&Regime=2]:3.0)[&Regime=0]:2.0, sp6[&Regime=3]:6.0)[&Regime=0]:1.0, sp7[&Regime=4]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStrNonUltraBig, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();

		// The values below are MLEs from mvMORPH
		double alphaFI = 7.142986;
		double alphaFM = 7.142986;
		double alphaRI = 7.84511;
		double alphaRM = 7.84511;
		
		// Setting weight matrices dimensions
		int n = myTree.getLeafNodeCount();
		ouWeight3FI = new double[n][6];
		ouWeight3FM = new double[n][5];
		ouWeight3RI = new double[n][6];
		ouWeight3RM = new double[n][5];
		
		OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 5, alphaFI, ouWeight3FI, false);
		OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 5, alphaFM, ouWeight3FM, true);
		OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 5, alphaRI, ouWeight3RI, false);
		OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 5, alphaRM, ouWeight3RM, true);
	}

	// Test 1
	@Test
	public void againstRweightMatTest1 () {
		
		Assert.assertEquals(7.815427e-28, 	ouWeight1FI[0][1], EPSILON);
		Assert.assertEquals(0, 				ouWeight1FI[2][2], EPSILON);
		Assert.assertEquals(1, 				ouWeight1FI[1][2], EPSILON);
		Assert.assertEquals(2.184887e-41, 	ouWeight1FI[1][0], EPSILON);
		
		Assert.assertEquals(1, 				ouWeight1FM[0][1], EPSILON);
		Assert.assertEquals(1, 				ouWeight1FM[2][2], EPSILON);
		Assert.assertEquals(0, 				ouWeight1FM[1][2], EPSILON);
		Assert.assertEquals(7.815427e-28, 	ouWeight1FM[1][0], EPSILON);
		
		Assert.assertEquals(2.208911e-38, 	ouWeight1RI[0][1], EPSILON);
		Assert.assertEquals(0, 				ouWeight1RI[2][2], EPSILON);
		Assert.assertEquals(1, 				ouWeight1RI[1][2], EPSILON);
		Assert.assertEquals(3.282974e-57, 	ouWeight1RI[1][0], EPSILON);
		
		Assert.assertEquals(1, 				ouWeight1RM[0][1], EPSILON);
		Assert.assertEquals(1, 				ouWeight1RM[2][2], EPSILON);
		Assert.assertEquals(0, 				ouWeight1RM[1][2], EPSILON);
		Assert.assertEquals(2.208911e-38, 	ouWeight1RM[1][0], EPSILON);
	}
	
	// Test 2
	@Test
	public void againstRweightMatTest2 () {
		   
		Assert.assertEquals(0.9495264, 	ouWeight2FI[0][1], EPSILON);
		Assert.assertEquals(0.8935126, 	ouWeight2FI[1][1], EPSILON);
		Assert.assertEquals(0.02392380, ouWeight2FI[2][0], EPSILON);
		Assert.assertEquals(0.10648738, ouWeight2FI[3][0], EPSILON);
		
		Assert.assertEquals(1, ouWeight2FM[0][0], EPSILON);
		Assert.assertEquals(1, ouWeight2FM[1][0], EPSILON);
		Assert.assertEquals(1, ouWeight2FM[2][0], EPSILON);
		Assert.assertEquals(1, ouWeight2FM[3][0], EPSILON);
		
		Assert.assertEquals(0.9680612, 	ouWeight2RI[0][1], EPSILON);
		Assert.assertEquals(0.9244492, 	ouWeight2RI[1][1], EPSILON);
		Assert.assertEquals(0.01350201, ouWeight2RI[2][0], EPSILON);
		Assert.assertEquals(0.07555081, ouWeight2RI[3][0], EPSILON);
		
		Assert.assertEquals(1, ouWeight2RM[0][0], EPSILON);
		Assert.assertEquals(1, ouWeight2RM[1][0], EPSILON);
		Assert.assertEquals(1, ouWeight2RM[2][0], EPSILON);
		Assert.assertEquals(1, ouWeight2RM[3][0], EPSILON);
	}
	
	// Test 3
	@Test
	public void againstRweightMatTest3 () {
		
		Assert.assertEquals(6.247136e-07, 	ouWeight3FI[0][1], EPSILON);
		Assert.assertEquals(0, 				ouWeight3FI[2][2], EPSILON);
		Assert.assertEquals(0.9999994, 		ouWeight3FI[1][2], EPSILON);
		Assert.assertEquals(1.927008e-22, 	ouWeight3FI[1][0], EPSILON);
		
		Assert.assertEquals(0.9999994, 		ouWeight3FM[0][1], EPSILON);
		Assert.assertEquals(0, 				ouWeight3FM[2][2], EPSILON);
		Assert.assertEquals(0, 				ouWeight3FM[1][2], EPSILON);
		Assert.assertEquals(6.247136e-07, 	ouWeight3FM[1][0], EPSILON);
		
		Assert.assertEquals(1.533995e-07, 	ouWeight3RI[0][1], EPSILON);
		Assert.assertEquals(0, 				ouWeight3RI[2][2], EPSILON);
		Assert.assertEquals(0.9999998, 		ouWeight3RI[1][2], EPSILON);
		Assert.assertEquals(1.413785e-24, 	ouWeight3RI[1][0], EPSILON);
		
		Assert.assertEquals(0.9999998, 		ouWeight3RM[0][1], EPSILON);
		Assert.assertEquals(0, 				ouWeight3RM[2][2], EPSILON);
		Assert.assertEquals(0, 				ouWeight3RM[1][2], EPSILON);
		Assert.assertEquals(1.533995e-07, 	ouWeight3RM[1][0], EPSILON);
	}
}
