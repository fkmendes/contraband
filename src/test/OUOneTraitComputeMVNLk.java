package test;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import beast.util.TreeParser;
import java.util.ArrayList;
import java.util.List;


import contraband.GeneralUtils;
import contraband.MVNUtils;
import beast.evolution.tree.Node;
import contraband.OUUtils;

public class OUOneTraitComputeMVNLk {
	
	final static double EPSILON = 1e-4;
	
	/* 
	 * The T matrix comes from the phylogenetic tree.
	 * The OU T matrix is equivalent to the V matrix in Butler & King (it's part of the variance 
	 * of the OU process).
	 * The W matrix is part of the mean of the OU process.
	 * 
	 * R: rootIsRandVar=true
	 * F: rootIsRandVar=false
	 * I: useRootMetaData=true
	 * M: useRootMetaData=false
	 */
		
	// We are using the same trees as in OUOneTraitcomputeOUTMatOneTraitTest.java)
	
	private static double ouLik1FI, ouLik1FM, ouLik1RI, ouLik1RM;
	private static double ouLik2FI, ouLik2FM, ouLik2RI, ouLik2RM;
	private static double ouLik3FI, ouLik3FM, ouLik3RI, ouLik3RM;
	
	// Test 1: Ultrametric tree (3 regimes)
	@BeforeClass
	public static void setUPTest1() throws Exception {

		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();
		
		double[] data = { 4.1, 4.5, 5.9, 0.0 };
		RealVector realData = new ArrayRealVector(data);
		
		// The values below are MLEs from mvMORPH
		double alphaFI = 31.20814;
		double alphaFM = 31.20814;
		double alphaRI = 43.35287;
		double alphaRM = 43.35287;
		
		double sigsqFI = 1.248328;		
		double sigsqFM = 1.248328;
		double sigsqRI = 1.734117;
		double sigsqRM = 1.734117;
		
		double[] thetaFI = {2.228585e-40, -4.047373e-16, 4.300000e+00, 5.900000e+00};
		double[] thetaFM = {8.128044e-27, 4.300000e+00, 5.900000e+00};
		double[] thetaRI = {3.348633e-56, -1.903330e-16, 4.300000e+00, 5.900000e+00};
		double[] thetaRM = {2.297268e-37, 4.300000e+00, 5.900000e+00};
		
		// Creating T matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
	
		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
	
		// OU T matrix (we multiply by sigma later, so not yet covariance matrix)
		double[][] ouTMatFI, ouTMatFM, ouTMatRI, ouTMatRM;
		ouTMatFI = new double[n][n];
		ouTMatFM = new double[n][n];
		ouTMatRI = new double[n][n];
		ouTMatRM = new double[n][n];

		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMatFI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMatFM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMatRI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMatRM, false);
		
		// Setting weight matrices dimensions
		double[][] ouWeightFI, ouWeightFM, ouWeightRI, ouWeightRM;
		ouWeightFI = new double[n][4];
		ouWeightFM = new double[n][3];
		ouWeightRI = new double[n][4];
		ouWeightRM = new double[n][3];
		 
		OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 3, alphaFI, ouWeightFI, false);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 3, alphaFM, ouWeightFM, true);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 3, alphaRI, ouWeightRI, false);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 3, alphaRM, ouWeightRM, true);
	    
	    // Preparing for likelihood computation
	    RealMatrix realCovFI = new Array2DRowRealMatrix(ouTMatFI);
		LUDecomposition realCovFILUD = new LUDecomposition(realCovFI);
		RealMatrix invRealCovFI = realCovFILUD.getSolver().getInverse(); 
		RealMatrix invFullCovFI = invRealCovFI.scalarMultiply(1/sigsqFI); 
		double varToNdetRealCovFI = Math.pow(sigsqFI, n) * realCovFILUD.getDeterminant(); 
		
	    RealMatrix realCovFM = new Array2DRowRealMatrix(ouTMatFM);
		LUDecomposition realCovFMLUD = new LUDecomposition(realCovFM);
		RealMatrix invRealCovFM = realCovFMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovFM = invRealCovFM.scalarMultiply(1/sigsqFM); 
		double varToNdetRealCovFM = Math.pow(sigsqFM, n) * realCovFMLUD.getDeterminant(); 
		
	    RealMatrix realCovRI = new Array2DRowRealMatrix(ouTMatRI);
		LUDecomposition realCovRILUD = new LUDecomposition(realCovRI);
		RealMatrix invRealCovRI = realCovRILUD.getSolver().getInverse(); 
		RealMatrix invFullCovRI = invRealCovRI.scalarMultiply(1/sigsqRI); 
		double varToNdetRealCovRI = Math.pow(sigsqRI, n) * realCovRILUD.getDeterminant(); 
		
	    RealMatrix realCovRM = new Array2DRowRealMatrix(ouTMatRM);
		LUDecomposition realCovRMLUD = new LUDecomposition(realCovRM);
		RealMatrix invRealCovRM = realCovRMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovRM = invRealCovRM.scalarMultiply(1/sigsqRM); 
		double varToNdetRealCovRM = Math.pow(sigsqRM, n) * realCovRMLUD.getDeterminant(); 
		
		RealVector realThetaFI = new ArrayRealVector(thetaFI);
		RealVector realThetaFM = new ArrayRealVector(thetaFM);
		RealVector realThetaRI = new ArrayRealVector(thetaRI);
		RealVector realThetaRM = new ArrayRealVector(thetaRM);	
		
		RealMatrix realOUWeightFI = new Array2DRowRealMatrix(ouWeightFI);
		RealMatrix realOUWeightFM = new Array2DRowRealMatrix(ouWeightFM);
		RealMatrix realOUWeightRI = new Array2DRowRealMatrix(ouWeightRI);
		RealMatrix realOUWeightRM = new Array2DRowRealMatrix(ouWeightRM);
		
//		System.out.println("sigsqFI " + sigsqFI);
//		System.out.println("realThetaFI: ");
//		GeneralUtils.displayRealVector(realOUWeightFI.operate(realThetaFI));
//		System.out.println("realData: ");
//		System.out.println(realData);
//		System.out.println("invFullCovFI");
//		GeneralUtils.displayRealMatrix(invFullCovFI);
//		System.out.println("varToNdetRealCovFI");
//		System.out.println(varToNdetRealCovFI);
//		System.out.println();
		
		// Finally computing likelihoods for OU
		ouLik1FI = MVNUtils.getMVNLk(n, sigsqFI, realOUWeightFI.operate(realThetaFI), realData, invFullCovFI, varToNdetRealCovFI);
		ouLik1FM = MVNUtils.getMVNLk(n, sigsqFM, realOUWeightFM.operate(realThetaFM), realData, invFullCovFM, varToNdetRealCovFM);
		ouLik1RI = MVNUtils.getMVNLk(n, sigsqRI, realOUWeightRI.operate(realThetaRI), realData, invFullCovRI, varToNdetRealCovRI);
		ouLik1RM = MVNUtils.getMVNLk(n, sigsqRM, realOUWeightRM.operate(realThetaRM), realData, invFullCovRM, varToNdetRealCovRM);
	}

	// Test 2: Non-ultrametric tree (1 regimes)
	@BeforeClass
	public static void setUPTest2() throws Exception {

		String treeStrNonUltra = "(((sp1[&Regime=0]:2.0, sp2[&Regime=0]:1.0)[&Regime=0]:1.0, sp3[&Regime=0]:4.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStrNonUltra, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();

		double[] data = {4.1, 4.5, 5.9, 0.0};
		RealVector realData = new ArrayRealVector(data);
		
		double alphaFI = 0.7465763;		// mvMORPH values
		double alphaFM = 1.40338e-08;
		double alphaRI = 0.8609833;
		double alphaRM = 0.7085376;
		
		double sigsqFI = 4.003551;	
		double sigsqFM = 1.237864;
		double sigsqRI = 4.601164;
		double sigsqRM = 6.841867;
		
		double[] thetaFI = {-33.591241, 6.449917};
		double[] thetaFM = {2.792045};
		double[] thetaRI = {-46.965464, 6.201598};
		double[] thetaRM = {3.586504};

		// Creating T matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
	
		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
	
		// OU T matrix (we multiply by sigma later, so not yet covariance matrix)
		double[][] ouTMatFI, ouTMatFM, ouTMatRI, ouTMatRM;
		ouTMatFI = new double[n][n];
		ouTMatFM = new double[n][n];
		ouTMatRI = new double[n][n];
		ouTMatRM = new double[n][n];

		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMatFI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMatFM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMatRI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMatRM, false);
	
		// Setting weight matrices dimensions
		double[][] ouWeightFI, ouWeightFM, ouWeightRI, ouWeightRM;
		ouWeightFI = new double[n][2]; // Don't forget to put the correct column dimensions (r = # regimes or r + 1 in the Isolated root case)
		ouWeightFM = new double[n][1];
		ouWeightRI = new double[n][2];
		ouWeightRM = new double[n][1];

	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 1, alphaFI, ouWeightFI, false);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 1, alphaFM, ouWeightFM, true);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 1, alphaRI, ouWeightRI, false);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 1, alphaRM, ouWeightRM, true);

	    // Preparing for likelihood computation
	    RealMatrix realCovFI = new Array2DRowRealMatrix(ouTMatFI);
		LUDecomposition realCovFILUD = new LUDecomposition(realCovFI);
		RealMatrix invRealCovFI = realCovFILUD.getSolver().getInverse(); 
		RealMatrix invFullCovFI = invRealCovFI.scalarMultiply(1/sigsqFI); 
		double varToNdetRealCovFI = Math.pow(sigsqFI, n) * realCovFILUD.getDeterminant(); 
		
	    RealMatrix realCovFM = new Array2DRowRealMatrix(ouTMatFM);
		LUDecomposition realCovFMLUD = new LUDecomposition(realCovFM);
		RealMatrix invRealCovFM = realCovFMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovFM = invRealCovFM.scalarMultiply(1/sigsqFM); 
		double varToNdetRealCovFM = Math.pow(sigsqFM, n) * realCovFMLUD.getDeterminant(); 
		
	    RealMatrix realCovRI = new Array2DRowRealMatrix(ouTMatRI);
		LUDecomposition realCovRILUD = new LUDecomposition(realCovRI);
		RealMatrix invRealCovRI = realCovRILUD.getSolver().getInverse(); 
		RealMatrix invFullCovRI = invRealCovRI.scalarMultiply(1/sigsqRI); 
		double varToNdetRealCovRI = Math.pow(sigsqRI, n) * realCovRILUD.getDeterminant(); 
		
	    RealMatrix realCovRM = new Array2DRowRealMatrix(ouTMatRM);
		LUDecomposition realCovRMLUD = new LUDecomposition(realCovRM);
		RealMatrix invRealCovRM = realCovRMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovRM = invRealCovRM.scalarMultiply(1/sigsqRM); 
		double varToNdetRealCovRM = Math.pow(sigsqRM, n) * realCovRMLUD.getDeterminant(); 
		
		RealVector realThetaFI = new ArrayRealVector(thetaFI);
		RealVector realThetaFM = new ArrayRealVector(thetaFM);
		RealVector realThetaRI = new ArrayRealVector(thetaRI);
		RealVector realThetaRM = new ArrayRealVector(thetaRM);
		
		RealMatrix realOUWeightFI = new Array2DRowRealMatrix(ouWeightFI);
		RealMatrix realOUWeightFM = new Array2DRowRealMatrix(ouWeightFM);
		RealMatrix realOUWeightRI = new Array2DRowRealMatrix(ouWeightRI);
		RealMatrix realOUWeightRM = new Array2DRowRealMatrix(ouWeightRM);

		// Finally computing likelihoods for OU
		ouLik2FI = MVNUtils.getMVNLk(n, sigsqFI, realOUWeightFI.operate(realThetaFI), realData, invFullCovFI, varToNdetRealCovFI);
		ouLik2FM = MVNUtils.getMVNLk(n, sigsqFM, realOUWeightFM.operate(realThetaFM), realData, invFullCovFM, varToNdetRealCovFM);
		ouLik2RI = MVNUtils.getMVNLk(n, sigsqRI, realOUWeightRI.operate(realThetaRI), realData, invFullCovRI, varToNdetRealCovRI);
		ouLik2RM = MVNUtils.getMVNLk(n, sigsqRM, realOUWeightRM.operate(realThetaRM), realData, invFullCovRM, varToNdetRealCovRM);
	}
	
	// Test 3: Non-ultrametric tree (5 regimes)
	@BeforeClass
	public static void setUPTest3() throws Exception {

		String treeStrNonUltraBig = "(((((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=0]:1.0)[&Regime=0]:2.0, (sp4[&Regime=2]:1.0, sp5[&Regime=2]:1.0)[&Regime=2]:3.0)[&Regime=0]:2.0, sp6[&Regime=3]:6.0)[&Regime=0]:1.0, sp7[&Regime=4]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStrNonUltraBig, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();
		
		double[] data = {4.1, 4.5, 5.9, 0.0, 3.2, 2.5, 5.4};
		RealVector realData = new ArrayRealVector(data);
		
		double alphaFI = 7.142986;		// mvMORPH values
		double alphaFM = 7.142986;
		double alphaRI = 7.84511;
		double alphaRM = 7.84511;
		
		double sigsqFI = 10.61289;	
		double sigsqFM = 10.61289;
		double sigsqRI = 11.65586;
		double sigsqRM = 11.65586;
		
		double[] thetaFI = {2.666338e-09, 5.900000e+00, 4.299999e+00, 1.600000e+00, 2.500000e+00, 5.400000e+00};
		double[] thetaFM = {5.900000e+00, 4.299999e+00, 1.600000e+00, 2.500000e+00, 5.400000e+00};
		double[] thetaRI = {3.244364e-10,  5.900000e+00, 4.300000e+00, 1.600000e+00, 2.500000e+00, 5.400000e+00};
		double[] thetaRM = {5.9, 4.3, 1.6, 2.5, 5.4};
		
		// Creating T matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
	
		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
	
		// OU T matrix (we multiply by sigma later, so not yet covariance matrix)
		double[][] ouTMatFI, ouTMatFM, ouTMatRI, ouTMatRM;
		ouTMatFI = new double[n][n];
		ouTMatFM = new double[n][n];
		ouTMatRI = new double[n][n];
		ouTMatRM = new double[n][n];

		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMatFI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMatFM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMatRI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMatRM, false);
	
		// Setting weight matrices dimensions
		double[][] ouWeightFI, ouWeightFM, ouWeightRI, ouWeightRM;
		ouWeightFI = new double[n][6]; // Don't forget to put the correct column dimensions
		ouWeightFM = new double[n][5];
		ouWeightRI = new double[n][6];
		ouWeightRM = new double[n][5];

	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 5, alphaFI, ouWeightFI, false);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 5, alphaFM, ouWeightFM, true);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 5, alphaRI, ouWeightRI, false);
	    OUUtils.computeWMatOneTrait(myTree, allLeafNodes, n, 5, alphaRM, ouWeightRM, true);
	    
	    // Preparing for likelihood computation
	    RealMatrix realCovFI = new Array2DRowRealMatrix(ouTMatFI);
		LUDecomposition realCovFILUD = new LUDecomposition(realCovFI);
		RealMatrix invRealCovFI = realCovFILUD.getSolver().getInverse(); 
		RealMatrix invFullCovFI = invRealCovFI.scalarMultiply(1/sigsqFI); 
		double varToNdetRealCovFI = Math.pow(sigsqFI, n) * realCovFILUD.getDeterminant(); 
		
	    RealMatrix realCovFM = new Array2DRowRealMatrix(ouTMatFM);
		LUDecomposition realCovFMLUD = new LUDecomposition(realCovFM);
		RealMatrix invRealCovFM = realCovFMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovFM = invRealCovFM.scalarMultiply(1/sigsqFM); 
		double varToNdetRealCovFM = Math.pow(sigsqFM, n) * realCovFMLUD.getDeterminant(); 
		
	    RealMatrix realCovRI = new Array2DRowRealMatrix(ouTMatRI);
		LUDecomposition realCovRILUD = new LUDecomposition(realCovRI);
		RealMatrix invRealCovRI = realCovRILUD.getSolver().getInverse(); 
		RealMatrix invFullCovRI = invRealCovRI.scalarMultiply(1/sigsqRI); 
		double varToNdetRealCovRI = Math.pow(sigsqRI, n) * realCovRILUD.getDeterminant(); 
		
	    RealMatrix realCovRM = new Array2DRowRealMatrix(ouTMatRM);
		LUDecomposition realCovRMLUD = new LUDecomposition(realCovRM);
		RealMatrix invRealCovRM = realCovRMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovRM = invRealCovRM.scalarMultiply(1/sigsqRM); 
		double varToNdetRealCovRM = Math.pow(sigsqRM, n) * realCovRMLUD.getDeterminant(); 
		
		RealVector realThetaFI = new ArrayRealVector(thetaFI);
		RealVector realThetaFM = new ArrayRealVector(thetaFM);
		RealVector realThetaRI = new ArrayRealVector(thetaRI);
		RealVector realThetaRM = new ArrayRealVector(thetaRM);	
		
		RealMatrix realOUWeightFI = new Array2DRowRealMatrix(ouWeightFI);
		RealMatrix realOUWeightFM = new Array2DRowRealMatrix(ouWeightFM);
		RealMatrix realOUWeightRI = new Array2DRowRealMatrix(ouWeightRI);
		RealMatrix realOUWeightRM = new Array2DRowRealMatrix(ouWeightRM);
		
//		System.out.println("Before");
//		
//		System.out.println("n = " + n);
//		System.out.println("sigsqFI = " + sigsqFI);
//		System.out.println("realThetaFI: ");
//		GeneralUtils.displayRealVector(realOUWeightFI.operate(realThetaFI));
//		System.out.println("realData: ");
//		System.out.println(realData);
//		System.out.println("invFullCovFI: ");
//		GeneralUtils.displayRealMatrix(invFullCovFI);
//		System.out.println("varToNdetRealCovFI");
//		System.out.println(varToNdetRealCovFI);
//		System.out.println();
		
		// Finally computing likelihoods for OU
		ouLik3FI = MVNUtils.getMVNLk(n, sigsqFI, realOUWeightFI.operate(realThetaFI), realData, invFullCovFI, varToNdetRealCovFI);	
		ouLik3FM = MVNUtils.getMVNLk(n, sigsqFM, realOUWeightFM.operate(realThetaFM), realData, invFullCovFM, varToNdetRealCovFM);
		ouLik3RI = MVNUtils.getMVNLk(n, sigsqRI, realOUWeightRI.operate(realThetaRI), realData, invFullCovRI, varToNdetRealCovRI);
		ouLik3RM = MVNUtils.getMVNLk(n, sigsqRM, realOUWeightRM.operate(realThetaRM), realData, invFullCovRM, varToNdetRealCovRM);

	}
	
	@Test // Test 1
	public void againstRcomputeLkTest1 () {
		Assert.assertEquals(2.148292, Math.log(ouLik1FI), EPSILON);
		Assert.assertEquals(2.148292, Math.log(ouLik1FM), EPSILON);
		Assert.assertEquals(2.148292, Math.log(ouLik1RI), EPSILON);
		Assert.assertEquals(2.148292, Math.log(ouLik1RM), EPSILON);	
	}
	
	@Test // Test 2
	public void againstRcomputeLkTest2 () {
		Assert.assertEquals(-7.630117, 	Math.log(ouLik2FI), EPSILON);
		Assert.assertEquals(-8.457486, 	Math.log(ouLik2FM), EPSILON);
		Assert.assertEquals(-7.63854, 	Math.log(ouLik2RI), EPSILON);
		Assert.assertEquals(-8.817273, 	Math.log(ouLik2RM), EPSILON);	
	}
	
	@Test // Test 3
	public void againstRcomputeLkTest3 () {
		Assert.assertEquals(-8.892192, 	Math.log(ouLik3FI), EPSILON);
		Assert.assertEquals(-8.892192, 	Math.log(ouLik3FM), EPSILON);
		Assert.assertEquals(-8.89219, 	Math.log(ouLik3RI), EPSILON);
		Assert.assertEquals(-8.89219, 	Math.log(ouLik3RM), EPSILON);
	}
	

}
