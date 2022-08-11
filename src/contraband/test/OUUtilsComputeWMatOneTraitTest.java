package contraband.test;

import java.util.List;

import beast.base.evolution.tree.Tree;
import contraband.clock.RateCategoryClockModel;
import contraband.clock.TreeToVCVMat;
import contraband.utils.GeneralUtils;
import contraband.utils.OUUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import contraband.clock.ColorManager;

/**
 * @author Pau Bravo and Fabio K. Mendes
 */

/*
 * This class contains unit tests for OUUtils
 */
public class OUUtilsComputeWMatOneTraitTest {

	final static double EPSILON = 1e-4;
	
	/* 
	 * The OU W mat is the "design" matrix in GLS lingo; it's part of the mean of the OU process
	 * (it gets multiplied by the vector of thetas in order to compute the mean).
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
	
	// We are using the same trees as in OUOneTraitcomputeOUTMatOneTraitTest.java)
	private static RealMatrix ouWeight1FI, ouWeight1FM, ouWeight1RI, ouWeight1RM;
	private static RealMatrix ouWeight2FI, ouWeight2FM, ouWeight2RI, ouWeight2RM;
	private static RealMatrix ouWeight3FI, ouWeight3FM, ouWeight3RI, ouWeight3RM;

	/*
	 * (1) Small ultrametric tree, in all combinations of R, F, I and M
	 */
	@Test
	public void testComputeOUWMatSmallTree() {
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();
		
		RealParameter colorValues = new RealParameter(new Double[] { 2.228585e-40, -4.047373e-16, 1.0, 1.0 }); // not used except for initialization of rcc (in OU, WMat is multiplied by this theta vector when computing lk)
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 1, 1, 2, 0, 1, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);
		// ColorManager colors = new ColorManager();
		// colors.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// The values below are MLEs from mvMORPH
		double alphaFI = 31.20814;
		double alphaFM = 31.20814;
		double alphaRI = 43.35287;
		double alphaRM = 43.35287;
		
		// Setting weight matrices dimensions
		int n = myTree.getLeafNodeCount();
		ouWeight1FI = new Array2DRowRealMatrix(n, 4);
		ouWeight1FM = new Array2DRowRealMatrix(n, 3);
		ouWeight1RI = new Array2DRowRealMatrix(n, 4);
		ouWeight1RM = new Array2DRowRealMatrix(n, 3);

		// filling WMat in place
	    OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 3, alphaFI, ouWeight1FI, true);
	    OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 3, alphaFM, ouWeight1FM, false);
	    OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 3, alphaRI, ouWeight1RI, true);
	    OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 3, alphaRM, ouWeight1RM, false);

		Assert.assertEquals(7.815427e-28, 	ouWeight1FI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(0, 			ouWeight1FI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(1, 			ouWeight1FI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(2.184887e-41, 	ouWeight1FI.getEntry(1, 0), EPSILON);

		Assert.assertEquals(1, 			ouWeight1FM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(1, 			ouWeight1FM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0, 			ouWeight1FM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(7.815427e-28, 	ouWeight1FM.getEntry(1, 0), EPSILON);

		Assert.assertEquals(2.208911e-38, 	ouWeight1RI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(0, 			ouWeight1RI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(1, 			ouWeight1RI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(3.282974e-57, 	ouWeight1RI.getEntry(1, 0), EPSILON);

		Assert.assertEquals(1, 			ouWeight1RM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(1, 			ouWeight1RM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0, 			ouWeight1RM.getEntry(1,	2), EPSILON);
		Assert.assertEquals(2.208911e-38, 	ouWeight1RM.getEntry(1, 0), EPSILON);
	}

	/*
	 * (2) Small non-ultrametric tree, in all combinations of R, F, I and M
	 */
	@Test
	public void testComputeOUWMatSmallTreeNonUltra() {
			
		String treeStr = "(((sp1:2.0,sp2:1.0):1.0,sp3:4.0):1.0,sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();

		RealParameter colorValues = new RealParameter(new Double[] { 0.1 }); // not used
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);
		// ColorManager colors = new ColorManager();
		// colors.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// The values below are MLEs from mvMORPH
		double alphaFI = 0.7465763;
		double alphaFM = 1.40338e-08;
		double alphaRI = 0.8609833;
		double alphaRM = 0.7085376;
			
		// Setting weight matrices dimensions
		int n = myTree.getLeafNodeCount();
		ouWeight2FI = new Array2DRowRealMatrix(n, 2);
		ouWeight2FM = new Array2DRowRealMatrix(n, 1);
		ouWeight2RI = new Array2DRowRealMatrix(n, 2);
		ouWeight2RM = new Array2DRowRealMatrix(n, 1);

		// filling WMat in place
		OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 1, alphaFI, ouWeight2FI, true);
		OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 1, alphaFM, ouWeight2FM, false);
		OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 1, alphaRI, ouWeight2RI, true);
		OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 1, alphaRM, ouWeight2RM, false);

		Assert.assertEquals(0.9495264, ouWeight2FI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(0.8935126, ouWeight2FI.getEntry(1, 1), EPSILON);
		Assert.assertEquals(0.02392380, ouWeight2FI.getEntry(2, 0), EPSILON);
		Assert.assertEquals(0.10648738, ouWeight2FI.getEntry(3, 0), EPSILON);

		Assert.assertEquals(1, ouWeight2FM.getEntry(0, 0), EPSILON);
		Assert.assertEquals(1, ouWeight2FM.getEntry(1, 0), EPSILON);
		Assert.assertEquals(1, ouWeight2FM.getEntry(2, 0), EPSILON);
		Assert.assertEquals(1, ouWeight2FM.getEntry(3, 0), EPSILON);

		Assert.assertEquals(0.9680612, ouWeight2RI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(0.9244492, ouWeight2RI.getEntry(1, 1), EPSILON);
		Assert.assertEquals(0.01350201, ouWeight2RI.getEntry(2, 0), EPSILON);
		Assert.assertEquals(0.07555081, ouWeight2RI.getEntry(3, 0), EPSILON);

		Assert.assertEquals(1, ouWeight2RM.getEntry(0, 0), EPSILON);
		Assert.assertEquals(1, ouWeight2RM.getEntry(1, 0), EPSILON);
		Assert.assertEquals(1, ouWeight2RM.getEntry(2, 0), EPSILON);
		Assert.assertEquals(1, ouWeight2RM.getEntry(3, 0), EPSILON);
	}

	/*
	 * (3) Larger non-ultrametric tree, in all combinations of R, F, I and M
	 */
	@Test
	public void testComputeOUWMatLargerTreeNonUltra() {

		String treeStr = "(((((sp1:1.0,sp2:1.0):1.0,sp3:1.0):2.0,(sp4:1.0,sp5:1.0):3.0):2.0,sp6:6.0):1.0,sp7:3.0);";
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();

		RealParameter colorValues = new RealParameter(new Double[] { 0.1, 0.1, 0.1, 0.1, 0.1 }); // not used
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 1, 1, 0, 2, 2, 3, 4, 1, 0, 2, 0, 1, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 5, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);
		// ColorManager colors = new ColorManager();
		// colors.initByName("nTraits", 1, "nColors", 5, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// The values below are MLEs from mvMORPH
		double alphaFI = 7.142986;
		double alphaFM = 7.142986;
		double alphaRI = 7.84511;
		double alphaRM = 7.84511;
			
		// Setting weight matrices dimensions
		int n = myTree.getLeafNodeCount();
		ouWeight3FI = new Array2DRowRealMatrix(n, 6);
		ouWeight3FM = new Array2DRowRealMatrix(n, 5);
		ouWeight3RI = new Array2DRowRealMatrix(n, 6);
		ouWeight3RM = new Array2DRowRealMatrix(n, 5);

		// filling WMat in place
		OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 5, alphaFI, ouWeight3FI, true);
		OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 5, alphaFM, ouWeight3FM, false);
		OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 5, alphaRI, ouWeight3RI, true);
		OUUtils.computeWMatOneTrait(colorAssignments.getValues(), myTree.getRoot(), allLeafNodes, n, 5, alphaRM, ouWeight3RM, false);

		Assert.assertEquals(6.247136e-07, 	ouWeight3FI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(0, 			ouWeight3FI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0.9999994, 	ouWeight3FI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.927008e-22, 	ouWeight3FI.getEntry(1, 0), EPSILON);

		Assert.assertEquals(0.9999994, 	ouWeight3FM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(0, 			ouWeight3FM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0, 			ouWeight3FM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(6.247136e-07, 	ouWeight3FM.getEntry(1, 0), EPSILON);

		Assert.assertEquals(1.533995e-07, 	ouWeight3RI.getEntry(0, 1), EPSILON);
		Assert.assertEquals(0, 			ouWeight3RI.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0.9999998, 	ouWeight3RI.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.413785e-24, 	ouWeight3RI.getEntry(1, 0), EPSILON);

		Assert.assertEquals(0.9999998, 	ouWeight3RM.getEntry(0, 1), EPSILON);
		Assert.assertEquals(0, 			ouWeight3RM.getEntry(2, 2), EPSILON);
		Assert.assertEquals(0, 			ouWeight3RM.getEntry(1, 2), EPSILON);
		Assert.assertEquals(1.533995e-07, 	ouWeight3RM.getEntry(1, 0), EPSILON);
	}
}
