package contraband.test;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Assert;
import org.junit.Test;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.Node;

import java.util.ArrayList;
import java.util.List;

import contraband.utils.MVNUtils;
import contraband.utils.OUUtils;

/**
 * @author Pau Bravo
 */

/*
 * This class contains unit tests for OUUtils
 */
public class OUMVNLikelihoodOneTraitFromUtilsTest {

	final static double EPSILON = 1e-4;

	/*
	 * The OU T matrix (after multiplying by sigma^2) is what is referred
	 * to in Butler & King as the V (variance) matrix of the OU process.
	 *
	 * We call it the T matrix here because initially we grab the evolutionary
	 * distances (without sigma^2), and then multiplied by sigma^2 in place --
	 * at this point, we are looking at the OU process V matrix.
	 * The OU W mat is the "design" matrix in GLS lingo; it's part of the mean
	 * of the OU process (it gets multiplied by the vector of thetas in order
	 * to compute the mean).
	 *
	 * The likelihood can be computed using the W matrix being tested here
	 * (as described in the Butler and King paper, say) but that's not how we do it
	 * in practice (in practice I traverse the tree, grabbing the theta and the length
	 * of each branch, and filling out the expectation of the OU process that way).
	 *
	 * This test and the implementation of the W matrix will remain in this package for
	 * legacy purposes (the first implementation used it).
	 *
	 * R: rootIsRandVar=true
	 * F: rootIsRandVar=false (assumes equilibrium distribution at root)
	 * I: useRootMetaData=true
	 * M: useRootMetaData=false (one fewer parameter)
	 */

	// We are using the same trees as in OUOneTraitcomputeOUTMatOneTraitTest.java
	private static double ouLik1FI, ouLik1FM, ouLik1RI, ouLik1RM;
	private static double ouLik2FI, ouLik2FM, ouLik2RI, ouLik2RM;
	private static double ouLik3FI, ouLik3FM, ouLik3RI, ouLik3RM;

	/*
	 * (1) Small ultrametric tree, in all combinations of R, F, I and M, with 3 optima
	 */
	@Test
	public void testComputeOULikSmallTreeThreeOpt() {

		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		Integer[] colorAssignments = new Integer[]{1, 1, 2, 0, 1, 0, 0, 0};
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();

		double[] data = {4.1, 4.5, 5.9, 0.0};
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

		double[] thetaFI = { 2.228585e-40, -4.047373e-16, 4.3, 5.9 };
		double[] thetaFM = { 8.128044e-27, 4.3, 5.9 };
		double[] thetaRI = { 3.348633e-56, -1.903330e-16, 4.3, 5.9 };
		double[] thetaRM = { 2.297268e-37, 4.3, 5.9 };

		// Creating T matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		String[] spNamesInPhyloTMatOrder = new String[n];

		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves, spNamesInPhyloTMatOrder);

		// OU T matrix (we multiply by sigma later, so not yet covariance matrix)
		RealMatrix ouTMatFI, ouTMatFM, ouTMatRI, ouTMatRM;
		ouTMatFI = new Array2DRowRealMatrix(n, n);
		ouTMatFM = new Array2DRowRealMatrix(n, n);
		ouTMatRI = new Array2DRowRealMatrix(n, n);
		ouTMatRM = new Array2DRowRealMatrix(n, n);

		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMatFI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMatFM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMatRI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMatRM, false);

		// Setting weight matrices dimensions
		RealMatrix ouWeightFI, ouWeightFM, ouWeightRI, ouWeightRM;
		ouWeightFI = new Array2DRowRealMatrix(n, 4);
		ouWeightFM = new Array2DRowRealMatrix(n, 3);
		ouWeightRI = new Array2DRowRealMatrix(n, 4);
		ouWeightRM = new Array2DRowRealMatrix(n, 3);

		OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 3, alphaFI, ouWeightFI, true);
		OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 3, alphaFM, ouWeightFM, false);
		OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 3, alphaRI, ouWeightRI, true);
		OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 3, alphaRM, ouWeightRM, false);

		// Preparing for likelihood computation

		LUDecomposition realCovFILUD = new LUDecomposition(ouTMatFI);
		RealMatrix invRealCovFI = realCovFILUD.getSolver().getInverse();
		RealMatrix invFullCovFI = invRealCovFI.scalarMultiply(1 / sigsqFI);
		double varToNdetRealCovFI = Math.pow(sigsqFI, n) * realCovFILUD.getDeterminant();

		LUDecomposition realCovFMLUD = new LUDecomposition(ouTMatFM);
		RealMatrix invRealCovFM = realCovFMLUD.getSolver().getInverse();
		RealMatrix invFullCovFM = invRealCovFM.scalarMultiply(1 / sigsqFM);
		double varToNdetRealCovFM = Math.pow(sigsqFM, n) * realCovFMLUD.getDeterminant();

		LUDecomposition realCovRILUD = new LUDecomposition(ouTMatRI);
		RealMatrix invRealCovRI = realCovRILUD.getSolver().getInverse();
		RealMatrix invFullCovRI = invRealCovRI.scalarMultiply(1 / sigsqRI);
		double varToNdetRealCovRI = Math.pow(sigsqRI, n) * realCovRILUD.getDeterminant();

		LUDecomposition realCovRMLUD = new LUDecomposition(ouTMatRM);
		RealMatrix invRealCovRM = realCovRMLUD.getSolver().getInverse();
		RealMatrix invFullCovRM = invRealCovRM.scalarMultiply(1 / sigsqRM);
		double varToNdetRealCovRM = Math.pow(sigsqRM, n) * realCovRMLUD.getDeterminant();

		RealVector realThetaFI = new ArrayRealVector(thetaFI);
		RealVector realThetaFM = new ArrayRealVector(thetaFM);
		RealVector realThetaRI = new ArrayRealVector(thetaRI);
		RealVector realThetaRM = new ArrayRealVector(thetaRM);

		// Finally computing likelihoods for OU
		ouLik1FI = MVNUtils.getMVNLk(n, ouWeightFI.operate(realThetaFI), realData, invFullCovFI, varToNdetRealCovFI);
		ouLik1FM = MVNUtils.getMVNLk(n, ouWeightFM.operate(realThetaFM), realData, invFullCovFM, varToNdetRealCovFM);
		ouLik1RI = MVNUtils.getMVNLk(n, ouWeightRI.operate(realThetaRI), realData, invFullCovRI, varToNdetRealCovRI);
		ouLik1RM = MVNUtils.getMVNLk(n, ouWeightRM.operate(realThetaRM), realData, invFullCovRM, varToNdetRealCovRM);

		Assert.assertEquals(2.148292, Math.log(ouLik1FI), EPSILON);
		Assert.assertEquals(2.148292, Math.log(ouLik1FM), EPSILON);
		Assert.assertEquals(2.148292, Math.log(ouLik1RI), EPSILON);
		Assert.assertEquals(2.148292, Math.log(ouLik1RM), EPSILON);
	}

	/*
	 * (2) Small non-ultrametric tree, in all combinations of R, F, I and M, with 1 optimum
	 */
	@Test
	public void testComputeOULikSmallTreeNonUltraOneOpt() {

		String treeStr = "(((sp1:2.0,sp2:1.0):1.0,sp3:4.0):1.0,sp4:3.0);";
		Integer[] colorAssignments = new Integer[]{ 0, 0, 0, 0, 0, 0, 0, 0 };
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();

		double[] data = { 4.1, 4.5, 5.9, 0.0 };
		RealVector realData = new ArrayRealVector(data);

		double alphaFI = 0.7465763;        // mvMORPH values
		double alphaFM = 1.40338e-08;
		double alphaRI = 0.8609833;
		double alphaRM = 0.7085376;

		double sigsqFI = 4.003551;
		double sigsqFM = 1.237864;
		double sigsqRI = 4.601164;
		double sigsqRM = 6.841867;

		double[] thetaFI = { -33.591241, 6.449917 };
		double[] thetaFM = { 2.792045};
		double[] thetaRI = { -46.965464, 6.201598 };
		double[] thetaRM = { 3.586504 };

		// Creating T matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		String[] spNamesInPhyloTMatOrder = new String[n];

		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves, spNamesInPhyloTMatOrder);

		// OU T matrix (we multiply by sigma later, so not yet covariance matrix)
		RealMatrix ouTMatFI, ouTMatFM, ouTMatRI, ouTMatRM;
		ouTMatFI = new Array2DRowRealMatrix(n, n);
		ouTMatFM = new Array2DRowRealMatrix(n, n);
		ouTMatRI = new Array2DRowRealMatrix(n, n);
		ouTMatRM = new Array2DRowRealMatrix(n, n);

		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMatFI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMatFM, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMatRI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMatRM, true);

		// Setting weight matrices dimensions
		RealMatrix ouWeightFI, ouWeightFM, ouWeightRI, ouWeightRM;
		ouWeightFI = new Array2DRowRealMatrix(n, 2); // Don't forget to put the correct column dimensions (r = # regimes or r + 1 in the Isolated root case)
		ouWeightFM = new Array2DRowRealMatrix(n, 1);
		ouWeightRI = new Array2DRowRealMatrix(n, 2);
		ouWeightRM = new Array2DRowRealMatrix(n, 1);

		OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 1, alphaFI, ouWeightFI, true);
		OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 1, alphaFM, ouWeightFM, false);
		OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 1, alphaRI, ouWeightRI, true);
		OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 1, alphaRM, ouWeightRM, false);

		// Preparing for likelihood computation
		LUDecomposition realCovFILUD = new LUDecomposition(ouTMatFI);
		RealMatrix invRealCovFI = realCovFILUD.getSolver().getInverse();
		RealMatrix invFullCovFI = invRealCovFI.scalarMultiply(1 / sigsqFI);
		double varToNdetRealCovFI = Math.pow(sigsqFI, n) * realCovFILUD.getDeterminant();

		LUDecomposition realCovFMLUD = new LUDecomposition(ouTMatFM);
		RealMatrix invRealCovFM = realCovFMLUD.getSolver().getInverse();
		RealMatrix invFullCovFM = invRealCovFM.scalarMultiply(1 / sigsqFM);
		double varToNdetRealCovFM = Math.pow(sigsqFM, n) * realCovFMLUD.getDeterminant();

		LUDecomposition realCovRILUD = new LUDecomposition(ouTMatRI);
		RealMatrix invRealCovRI = realCovRILUD.getSolver().getInverse();
		RealMatrix invFullCovRI = invRealCovRI.scalarMultiply(1 / sigsqRI);
		double varToNdetRealCovRI = Math.pow(sigsqRI, n) * realCovRILUD.getDeterminant();

		LUDecomposition realCovRMLUD = new LUDecomposition(ouTMatRM);
		RealMatrix invRealCovRM = realCovRMLUD.getSolver().getInverse();
		RealMatrix invFullCovRM = invRealCovRM.scalarMultiply(1 / sigsqRM);
		double varToNdetRealCovRM = Math.pow(sigsqRM, n) * realCovRMLUD.getDeterminant();

		RealVector realThetaFI = new ArrayRealVector(thetaFI);
		RealVector realThetaFM = new ArrayRealVector(thetaFM);
		RealVector realThetaRI = new ArrayRealVector(thetaRI);
		RealVector realThetaRM = new ArrayRealVector(thetaRM);

		// Finally computing likelihoods for OU
		ouLik2FI = MVNUtils.getMVNLk(n, ouWeightFI.operate(realThetaFI), realData, invFullCovFI, varToNdetRealCovFI);
		ouLik2FM = MVNUtils.getMVNLk(n, ouWeightFM.operate(realThetaFM), realData, invFullCovFM, varToNdetRealCovFM);
		ouLik2RI = MVNUtils.getMVNLk(n, ouWeightRI.operate(realThetaRI), realData, invFullCovRI, varToNdetRealCovRI);
		ouLik2RM = MVNUtils.getMVNLk(n, ouWeightRM.operate(realThetaRM), realData, invFullCovRM, varToNdetRealCovRM);

		Assert.assertEquals(-7.630117, 	Math.log(ouLik2FI), EPSILON);
		Assert.assertEquals(-8.457486, 	Math.log(ouLik2FM), EPSILON);
		Assert.assertEquals(-7.63854, 	Math.log(ouLik2RI), EPSILON);
		Assert.assertEquals(-8.817273, 	Math.log(ouLik2RM), EPSILON);
	}

	/*
	 * (3) Large non-ultrametric tree, in all combinations of R, F, I and M, with 1 optimum
	 */
	@Test
	public void testComputeOULikLargeTreeNonUltraFiveOpt() {

		String treeStr = "(((((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=0]:1.0)[&Regime=0]:2.0, (sp4[&Regime=2]:1.0, sp5[&Regime=2]:1.0)[&Regime=2]:3.0)[&Regime=0]:2.0, sp6[&Regime=3]:6.0)[&Regime=0]:1.0, sp7[&Regime=4]:3.0)[&Regime=0];";
		Integer[] colorAssignments = new Integer[] { 1, 1, 0, 2, 2, 3, 4, 0, 1, 0, 2, 0, 0, 0 };
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();

		double[] data = { 4.1, 4.5, 5.9, 0.0, 3.2, 2.5, 5.4 };
		RealVector realData = new ArrayRealVector(data);

		double alphaFI = 7.142986;		// mvMORPH values
		double alphaFM = 7.142986;
		double alphaRI = 7.84511;
		double alphaRM = 7.84511;

		double sigsqFI = 10.61289;
		double sigsqFM = 10.61289;
		double sigsqRI = 11.65586;
		double sigsqRM = 11.65586;

		double[] thetaFI = { 2.666338e-09, 5.9, 4.299999, 1.6, 2.500000e+00, 5.4 };
		double[] thetaFM = { 5.9, 4.299999, 1.6, 2.5, 5.4 };
		double[] thetaRI = { 3.244364e-10,  5.9, 4.3, 1.6, 2.5, 5.4 };
		double[] thetaRM = { 5.9, 4.3, 1.6, 2.5, 5.4 };

		// Creating T matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		String[] spNamesInPhyloTMatOrder = new String[n];

		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves, spNamesInPhyloTMatOrder);

		// OU T matrix (we multiply by sigma later, so not yet covariance matrix)
		RealMatrix ouTMatFI, ouTMatFM, ouTMatRI, ouTMatRM;
		ouTMatFI = new Array2DRowRealMatrix(n, n);
		ouTMatFM = new Array2DRowRealMatrix(n, n);
		ouTMatRI = new Array2DRowRealMatrix(n, n);
		ouTMatRM = new Array2DRowRealMatrix(n, n);

		OUUtils.computeOUTMatOneTrait(n, alphaFI, tMatInput, ouTMatFI, true);
		OUUtils.computeOUTMatOneTrait(n, alphaFM, tMatInput, ouTMatFM, true);
		OUUtils.computeOUTMatOneTrait(n, alphaRI, tMatInput, ouTMatRI, false);
		OUUtils.computeOUTMatOneTrait(n, alphaRM, tMatInput, ouTMatRM, false);

		// Setting weight matrices dimensions
		RealMatrix ouWeightFI, ouWeightFM, ouWeightRI, ouWeightRM;
		ouWeightFI = new Array2DRowRealMatrix(n, 6); // Don't forget to put the correct column dimensions
		ouWeightFM = new Array2DRowRealMatrix(n, 5);
		ouWeightRI = new Array2DRowRealMatrix(n, 6);
		ouWeightRM = new Array2DRowRealMatrix(n, 5);

	    OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 5, alphaFI, ouWeightFI, true);
	    OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 5, alphaFM, ouWeightFM, false);
	    OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 5, alphaRI, ouWeightRI, true);
	    OUUtils.computeWMatOneTrait(colorAssignments, myTree.getRoot(), allLeafNodes, n, 5, alphaRM, ouWeightRM, false);

	    // Preparing for likelihood computation
	    LUDecomposition realCovFILUD = new LUDecomposition(ouTMatFI);
		RealMatrix invRealCovFI = realCovFILUD.getSolver().getInverse();
		RealMatrix invFullCovFI = invRealCovFI.scalarMultiply(1/sigsqFI);
		double varToNdetRealCovFI = Math.pow(sigsqFI, n) * realCovFILUD.getDeterminant();

		LUDecomposition realCovFMLUD = new LUDecomposition(ouTMatFM);
		RealMatrix invRealCovFM = realCovFMLUD.getSolver().getInverse();
		RealMatrix invFullCovFM = invRealCovFM.scalarMultiply(1/sigsqFM);
		double varToNdetRealCovFM = Math.pow(sigsqFM, n) * realCovFMLUD.getDeterminant();

		LUDecomposition realCovRILUD = new LUDecomposition(ouTMatRI);
		RealMatrix invRealCovRI = realCovRILUD.getSolver().getInverse();
		RealMatrix invFullCovRI = invRealCovRI.scalarMultiply(1/sigsqRI);
		double varToNdetRealCovRI = Math.pow(sigsqRI, n) * realCovRILUD.getDeterminant();

		LUDecomposition realCovRMLUD = new LUDecomposition(ouTMatRM);
		RealMatrix invRealCovRM = realCovRMLUD.getSolver().getInverse();
		RealMatrix invFullCovRM = invRealCovRM.scalarMultiply(1/sigsqRM);
		double varToNdetRealCovRM = Math.pow(sigsqRM, n) * realCovRMLUD.getDeterminant();

		RealVector realThetaFI = new ArrayRealVector(thetaFI);
		RealVector realThetaFM = new ArrayRealVector(thetaFM);
		RealVector realThetaRI = new ArrayRealVector(thetaRI);
		RealVector realThetaRM = new ArrayRealVector(thetaRM);

		// Finally computing likelihoods for OU
		ouLik3FI = MVNUtils.getMVNLk(n, ouWeightFI.operate(realThetaFI), realData, invFullCovFI, varToNdetRealCovFI);
		ouLik3FM = MVNUtils.getMVNLk(n, ouWeightFM.operate(realThetaFM), realData, invFullCovFM, varToNdetRealCovFM);
		ouLik3RI = MVNUtils.getMVNLk(n, ouWeightRI.operate(realThetaRI), realData, invFullCovRI, varToNdetRealCovRI);
		ouLik3RM = MVNUtils.getMVNLk(n, ouWeightRM.operate(realThetaRM), realData, invFullCovRM, varToNdetRealCovRM);

		Assert.assertEquals(-8.892192, 	Math.log(ouLik3FI), EPSILON);
		Assert.assertEquals(-8.892192, 	Math.log(ouLik3FM), EPSILON);
		Assert.assertEquals(-8.89219, 	Math.log(ouLik3RI), EPSILON);
		Assert.assertEquals(-8.89219, 	Math.log(ouLik3RM), EPSILON);
	}
}