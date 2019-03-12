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

public class OUOneTraitcomputeMVNLk {
	
	final static double EPSILON = 1e-4;
	
	//For the covariance matrices we can assume two different hipothesis:
	//'F' suffix refers to assuming a fixed parameter Root
	//'R' suffix refers to assuming a random variable Root (stationary distribution)
	//For the weight matrices we can assume two different hipothesis:
	//'I' suffix refers to isolating the root optimum value in the weight Matrix
	//'M' suffix refers to merging the root parameter with the optimum parameter associated with the eldest selective regime

	//Every test has four different outputs according to the combination of the previous situations: FI, FM, RI, RM
		
	// 1, 2, 3 refers to Test 1, Test 2 and Test 3 respectively (WE ARE USING THE SAME TREES AS IN OUOneTraitcomputeOUTMatOneTraitTest.java)
	
	private static double ouLik1FI, ouLik1FM, ouLik1RI, ouLik1RM;
	private static double ouLik2FI, ouLik2FM, ouLik2RI, ouLik2RM;
	private static double ouLik3FI, ouLik3FM, ouLik3RI, ouLik3RM;
	
	@BeforeClass	// Test 1: Ultrametric tree (3 regimes)
	public static void setUPTest1() throws Exception {

		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		double[] data = {4.1, 4.5, 5.9, 0.0};
		
		double alphaFI = 31.20814;				// mvMORPH values
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
		
			// Creating phylogenetic covariance matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
	
		MVNUtils.populateVcvMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
	
			// OU covariance matrices
		double[][] ouCovFI, ouCovFM, ouCovRI, ouCovRM;
		ouCovFI = new double[n][n];
		ouCovFM = new double[n][n];
		ouCovRI = new double[n][n];
		ouCovRM = new double[n][n];

		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouCovFI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouCovFM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouCovRI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouCovRM, false);
	
		// We DON't multiply the covariance matrices by sigma because in the likelihood calculation we want
		// want to hava them separately for higher performance purposes
		
			// Setting weight matrices dimensions
		double[][] ouWeightFI, ouWeightFM, ouWeightRI, ouWeightRM;
		ouWeightFI = new double[n][4];
		ouWeightFM = new double[n][3];
		ouWeightRI = new double[n][4];
		ouWeightRM = new double[n][3];
			// computing the values for all cases
	    OUUtils.computeWMatOneTrait(myTree, n, 3, alphaFI, ouWeightFI, false);
	    OUUtils.computeWMatOneTrait(myTree, n, 3, alphaFM, ouWeightFM, true);
	    OUUtils.computeWMatOneTrait(myTree, n, 3, alphaRI, ouWeightRI, false);
	    OUUtils.computeWMatOneTrait(myTree, n, 3, alphaRM, ouWeightRM, true);
	    
//		GeneralUtils.display2DArray(ouWeightFI);
//		System.out.println();
	    
	    	// Likelihood computation
	    RealMatrix realCovFI = new Array2DRowRealMatrix(ouCovFI);
		LUDecomposition realCovFILUD = new LUDecomposition(realCovFI);
		RealMatrix invRealCovFI = realCovFILUD.getSolver().getInverse(); 
		RealMatrix invFullCovFI = invRealCovFI.scalarMultiply(1/sigsqFI); 
		double varToNdetRealCovFI = Math.pow(sigsqFI, n) * realCovFILUD.getDeterminant(); 
		
	    RealMatrix realCovFM = new Array2DRowRealMatrix(ouCovFM);
		LUDecomposition realCovFMLUD = new LUDecomposition(realCovFM);
		RealMatrix invRealCovFM = realCovFMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovFM = invRealCovFM.scalarMultiply(1/sigsqFM); 
		double varToNdetRealCovFM = Math.pow(sigsqFM, n) * realCovFMLUD.getDeterminant(); 
		
	    RealMatrix realCovRI = new Array2DRowRealMatrix(ouCovRI);
		LUDecomposition realCovRILUD = new LUDecomposition(realCovRI);
		RealMatrix invRealCovRI = realCovRILUD.getSolver().getInverse(); 
		RealMatrix invFullCovRI = invRealCovRI.scalarMultiply(1/sigsqRI); 
		double varToNdetRealCovRI = Math.pow(sigsqRI, n) * realCovRILUD.getDeterminant(); 
		
	    RealMatrix realCovRM = new Array2DRowRealMatrix(ouCovRM);
		LUDecomposition realCovRMLUD = new LUDecomposition(realCovRM);
		RealMatrix invRealCovRM = realCovRMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovRM = invRealCovRM.scalarMultiply(1/sigsqRM); 
		double varToNdetRealCovRM = Math.pow(sigsqRM, n) * realCovRMLUD.getDeterminant(); 
		
		RealVector realData = new ArrayRealVector(data);
		
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
		
		ouLik1FI = MVNUtils.computeMVNLk(n, sigsqFI, realOUWeightFI.operate(realThetaFI), realData, invFullCovFI, varToNdetRealCovFI);
		ouLik1FM = MVNUtils.computeMVNLk(n, sigsqFM, realOUWeightFM.operate(realThetaFM), realData, invFullCovFM, varToNdetRealCovFM);
		ouLik1RI = MVNUtils.computeMVNLk(n, sigsqRI, realOUWeightRI.operate(realThetaRI), realData, invFullCovRI, varToNdetRealCovRI);
		ouLik1RM = MVNUtils.computeMVNLk(n, sigsqRM, realOUWeightRM.operate(realThetaRM), realData, invFullCovRM, varToNdetRealCovRM);
	}

	@BeforeClass	// Test 2: non Ultrametric tree (1 regimes)
	public static void setUPTest2() throws Exception {

		String treeStrNonUltra = "(((sp1[&Regime=0]:2.0, sp2[&Regime=0]:1.0)[&Regime=0]:1.0, sp3[&Regime=0]:4.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStrNonUltra, false, false, true, 0);
		
		double[] data = {4.1, 4.5, 5.9, 0.0};
		
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

			// Creating phylogenetic covariance matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
	
		MVNUtils.populateVcvMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
	
			// OU covariance matrices
		double[][] ouCovFI, ouCovFM, ouCovRI, ouCovRM;
		ouCovFI = new double[n][n];
		ouCovFM = new double[n][n];
		ouCovRI = new double[n][n];
		ouCovRM = new double[n][n];

		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouCovFI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouCovFM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouCovRI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouCovRM, false);
	
		// We DON't multiply the covariance matrices by sigma because in the likelihood calculation we want
		// want to hava them separately for higher performance purposes


			// Setting weight matrices dimensions
		double[][] ouWeightFI, ouWeightFM, ouWeightRI, ouWeightRM;
		ouWeightFI = new double[n][2]; // Don't forget to put the correct column dimensions (r = # regimes or r + 1 in the Isolated root case)
		ouWeightFM = new double[n][1];
		ouWeightRI = new double[n][2];
		ouWeightRM = new double[n][1];
			// computing the values for all cases
	    OUUtils.computeWMatOneTrait(myTree, n, 1, alphaFI, ouWeightFI, false);
	    OUUtils.computeWMatOneTrait(myTree, n, 1, alphaFM, ouWeightFM, true);
	    OUUtils.computeWMatOneTrait(myTree, n, 1, alphaRI, ouWeightRI, false);
	    OUUtils.computeWMatOneTrait(myTree, n, 1, alphaRM, ouWeightRM, true);

	    	// Likelihood computation
	    RealMatrix realCovFI = new Array2DRowRealMatrix(ouCovFI);
		LUDecomposition realCovFILUD = new LUDecomposition(realCovFI);
		RealMatrix invRealCovFI = realCovFILUD.getSolver().getInverse(); 
		RealMatrix invFullCovFI = invRealCovFI.scalarMultiply(1/sigsqFI); 
		double varToNdetRealCovFI = Math.pow(sigsqFI, n) * realCovFILUD.getDeterminant(); 
		
	    RealMatrix realCovFM = new Array2DRowRealMatrix(ouCovFM);
		LUDecomposition realCovFMLUD = new LUDecomposition(realCovFM);
		RealMatrix invRealCovFM = realCovFMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovFM = invRealCovFM.scalarMultiply(1/sigsqFM); 
		double varToNdetRealCovFM = Math.pow(sigsqFM, n) * realCovFMLUD.getDeterminant(); 
		
	    RealMatrix realCovRI = new Array2DRowRealMatrix(ouCovRI);
		LUDecomposition realCovRILUD = new LUDecomposition(realCovRI);
		RealMatrix invRealCovRI = realCovRILUD.getSolver().getInverse(); 
		RealMatrix invFullCovRI = invRealCovRI.scalarMultiply(1/sigsqRI); 
		double varToNdetRealCovRI = Math.pow(sigsqRI, n) * realCovRILUD.getDeterminant(); 
		
	    RealMatrix realCovRM = new Array2DRowRealMatrix(ouCovRM);
		LUDecomposition realCovRMLUD = new LUDecomposition(realCovRM);
		RealMatrix invRealCovRM = realCovRMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovRM = invRealCovRM.scalarMultiply(1/sigsqRM); 
		double varToNdetRealCovRM = Math.pow(sigsqRM, n) * realCovRMLUD.getDeterminant(); 

		RealVector realData = new ArrayRealVector(data);
		
		RealVector realThetaFI = new ArrayRealVector(thetaFI);
		RealVector realThetaFM = new ArrayRealVector(thetaFM);
		RealVector realThetaRI = new ArrayRealVector(thetaRI);
		RealVector realThetaRM = new ArrayRealVector(thetaRM);
		
		RealMatrix realOUWeightFI = new Array2DRowRealMatrix(ouWeightFI);
		RealMatrix realOUWeightFM = new Array2DRowRealMatrix(ouWeightFM);
		RealMatrix realOUWeightRI = new Array2DRowRealMatrix(ouWeightRI);
		RealMatrix realOUWeightRM = new Array2DRowRealMatrix(ouWeightRM);

		ouLik2FI = MVNUtils.computeMVNLk(n, sigsqFI, realOUWeightFI.operate(realThetaFI), realData, invFullCovFI, varToNdetRealCovFI);
		ouLik2FM = MVNUtils.computeMVNLk(n, sigsqFM, realOUWeightFM.operate(realThetaFM), realData, invFullCovFM, varToNdetRealCovFM);
		ouLik2RI = MVNUtils.computeMVNLk(n, sigsqRI, realOUWeightRI.operate(realThetaRI), realData, invFullCovRI, varToNdetRealCovRI);
		ouLik2RM = MVNUtils.computeMVNLk(n, sigsqRM, realOUWeightRM.operate(realThetaRM), realData, invFullCovRM, varToNdetRealCovRM);
	}
	
	@BeforeClass	// Test 3: non Ultrametric tree (5 regimes)
	public static void setUPTest3() throws Exception {

		String treeStrNonUltraBig = "(((((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=0]:1.0)[&Regime=0]:2.0, (sp4[&Regime=2]:1.0, sp5[&Regime=2]:1.0)[&Regime=2]:3.0)[&Regime=0]:2.0, sp6[&Regime=3]:6.0)[&Regime=0]:1.0, sp7[&Regime=4]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStrNonUltraBig, false, false, true, 0);
		
		double[] data = {4.1, 4.5, 5.9, 0.0, 3.2, 2.5, 5.4};
		
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
		
			// Creating phylogenetic covariance matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
	
		MVNUtils.populateVcvMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
	
			// OU covariance matrices
		double[][] ouCovFI, ouCovFM, ouCovRI, ouCovRM;
		ouCovFI = new double[n][n];
		ouCovFM = new double[n][n];
		ouCovRI = new double[n][n];
		ouCovRM = new double[n][n];

		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouCovFI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouCovFM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouCovRI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouCovRM, false);
	
		// We DON't multiply the covariance matrices by sigma because in the likelihood calculation we want
		// want to have them separately for higher performance purposes
		

			// Setting weight matrices dimensions
		double[][] ouWeightFI, ouWeightFM, ouWeightRI, ouWeightRM;
		ouWeightFI = new double[n][6]; // Don't forget to put the correct column dimensions
		ouWeightFM = new double[n][5];
		ouWeightRI = new double[n][6];
		ouWeightRM = new double[n][5];
			// computing the values for all cases
	    OUUtils.computeWMatOneTrait(myTree, n, 5, alphaFI, ouWeightFI, false);
	    OUUtils.computeWMatOneTrait(myTree, n, 5, alphaFM, ouWeightFM, true);
	    OUUtils.computeWMatOneTrait(myTree, n, 5, alphaRI, ouWeightRI, false);
	    OUUtils.computeWMatOneTrait(myTree, n, 5, alphaRM, ouWeightRM, true);
	    
	    	// Likelihood computation
	    RealMatrix realCovFI = new Array2DRowRealMatrix(ouCovFI);
		LUDecomposition realCovFILUD = new LUDecomposition(realCovFI);
		RealMatrix invRealCovFI = realCovFILUD.getSolver().getInverse(); 
		RealMatrix invFullCovFI = invRealCovFI.scalarMultiply(1/sigsqFI); 
		double varToNdetRealCovFI = Math.pow(sigsqFI, n) * realCovFILUD.getDeterminant(); 
		
	    RealMatrix realCovFM = new Array2DRowRealMatrix(ouCovFM);
		LUDecomposition realCovFMLUD = new LUDecomposition(realCovFM);
		RealMatrix invRealCovFM = realCovFMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovFM = invRealCovFM.scalarMultiply(1/sigsqFM); 
		double varToNdetRealCovFM = Math.pow(sigsqFM, n) * realCovFMLUD.getDeterminant(); 
		
	    RealMatrix realCovRI = new Array2DRowRealMatrix(ouCovRI);
		LUDecomposition realCovRILUD = new LUDecomposition(realCovRI);
		RealMatrix invRealCovRI = realCovRILUD.getSolver().getInverse(); 
		RealMatrix invFullCovRI = invRealCovRI.scalarMultiply(1/sigsqRI); 
		double varToNdetRealCovRI = Math.pow(sigsqRI, n) * realCovRILUD.getDeterminant(); 
		
	    RealMatrix realCovRM = new Array2DRowRealMatrix(ouCovRM);
		LUDecomposition realCovRMLUD = new LUDecomposition(realCovRM);
		RealMatrix invRealCovRM = realCovRMLUD.getSolver().getInverse(); 
		RealMatrix invFullCovRM = invRealCovRM.scalarMultiply(1/sigsqRM); 
		double varToNdetRealCovRM = Math.pow(sigsqRM, n) * realCovRMLUD.getDeterminant(); 
		
		RealVector realData = new ArrayRealVector(data);
		
		RealVector realThetaFI = new ArrayRealVector(thetaFI);
		RealVector realThetaFM = new ArrayRealVector(thetaFM);
		RealVector realThetaRI = new ArrayRealVector(thetaRI);
		RealVector realThetaRM = new ArrayRealVector(thetaRM);	
		
		RealMatrix realOUWeightFI = new Array2DRowRealMatrix(ouWeightFI);
		RealMatrix realOUWeightFM = new Array2DRowRealMatrix(ouWeightFM);
		RealMatrix realOUWeightRI = new Array2DRowRealMatrix(ouWeightRI);
		RealMatrix realOUWeightRM = new Array2DRowRealMatrix(ouWeightRM);
		
		System.out.println("Before");
		
		System.out.println("n = " + n);
		System.out.println("sigsqFI = " + sigsqFI);
		System.out.println("realThetaFI: ");
		GeneralUtils.displayRealVector(realOUWeightFI.operate(realThetaFI));
		System.out.println("realData: ");
		System.out.println(realData);
		System.out.println("invFullCovFI: ");
		GeneralUtils.displayRealMatrix(invFullCovFI);
		System.out.println("varToNdetRealCovFI");
		System.out.println(varToNdetRealCovFI);
		System.out.println();
		

		ouLik3FI = MVNUtils.computeMVNLk(n, sigsqFI, realOUWeightFI.operate(realThetaFI), realData, invFullCovFI, varToNdetRealCovFI);
		
		ouLik3FM = MVNUtils.computeMVNLk(n, sigsqFM, realOUWeightFM.operate(realThetaFM), realData, invFullCovFM, varToNdetRealCovFM);
		ouLik3RI = MVNUtils.computeMVNLk(n, sigsqRI, realOUWeightRI.operate(realThetaRI), realData, invFullCovRI, varToNdetRealCovRI);
		ouLik3RM = MVNUtils.computeMVNLk(n, sigsqRM, realOUWeightRM.operate(realThetaRM), realData, invFullCovRM, varToNdetRealCovRM);

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
