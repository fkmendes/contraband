package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.ColorManager;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;
import contraband.RateCategoryClockModel;
import contraband.TreeToVCVMat;

/*
 * Small tree, three optima, with and without stationary distr for root value; having root value as a separate
 * parameter or setting it to be optimum with index 0 (the index 0, and the optima have to be in increasing order,
 * this can have effects when comparing mvMORPH's results and ours, so make sure they match)
 */
public class OUMVNLikelihoodOneTraitTest1 {

	double lnLk1, lnLk2, lnLk3, lnLk4;
	final static double EPSILON = 1e-6;
	
	@Before
	public void setUp() throws Exception {
		// tree
		// String treeStr = "(((sp1[&Regime=1]:1.0,sp2[&Regime=1]:1.0)[&Regime=1]:1.0,sp3[&Regime=2]:2.0)[&Regime=0]:1.0,sp4[&Regime=0]:3.0)[&Regime=0];";
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data		
		// RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9, 0.0 });
		RealParameter oneTraitValues = new RealParameter(new Double[] { 0.237649365136715, 0.295018750722361, 0.881225138279161, 0.206222932069516 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
		
		// thetas
		RealParameter colorValues1 = new RealParameter(new Double[] { 0.206222932117995, 0.26633408087427, 0.88122539543514 });
		RealParameter colorValues2 = new RealParameter(new Double[] { 0.206222932069516, 0.266334080825641, 0.881225395384966 });
		RealParameter colorValues3 = new RealParameter(new Double[] { 0.206222932069532, 0.266334058036916, 0.881225139484772 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 1, 1, 2, 0, 1, 0, 0 });
				
		RateCategoryClockModel rcc1 = new RateCategoryClockModel();
		rcc1.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues1, "tree", myTree);
		RateCategoryClockModel rcc2 = new RateCategoryClockModel();
		rcc2.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues2, "tree", myTree);
		RateCategoryClockModel rcc3 = new RateCategoryClockModel();
		rcc3.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues3, "tree", myTree);
		
		TreeToVCVMat optima1 = new TreeToVCVMat();
		optima1.initByName("branchRateModel", rcc1, "tree", myTree, "coalCorrection", false);
		TreeToVCVMat optima2 = new TreeToVCVMat();
		optima2.initByName("branchRateModel", rcc2, "tree", myTree, "coalCorrection", false);
		TreeToVCVMat optima3 = new TreeToVCVMat();
		optima3.initByName("branchRateModel", rcc3, "tree", myTree, "coalCorrection", false);
		TreeToVCVMat optima4 = new TreeToVCVMat();
		optima4.initByName("branchRateModel", rcc3, "tree", myTree, "coalCorrection", false);
		
		/*
		 * Old parameterization using W matrix and ColorManager
		 */
//		ColorManager optima1 = new ColorManager();
//		optima1.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues1, "colorAssignments", colorAssignments1, "coalCorrection", false);
//		ColorManager optima2 = new ColorManager();
//		optima2.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues2, "colorAssignments", colorAssignments1, "coalCorrection", false);
//		ColorManager optima3 = new ColorManager();
//		optima3.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues3, "colorAssignments", colorAssignments1, "coalCorrection", false);
//		ColorManager optima4 = new ColorManager();
//		optima4.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues3, "colorAssignments", colorAssignments1, "coalCorrection", false);
		
		// sigmasq
		// Double[] sigmasqInput = new Double[] { 1.248328 };
		Double[] sigmasqInput = new Double[] { 0.006082604 };
		RealParameter sigmasq1 = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 0.008287661 };
		RealParameter sigmasq2 = new RealParameter(sigmasqInput);
		
		// alpha
		// Double[] alphaInput = new Double[] { 31.20814 };
		Double[] alphaInput = new Double[] { 7.390366 };
		RealParameter alpha1 = new RealParameter(alphaInput);
		alphaInput = new Double[] { 10.07163 };
		RealParameter alpha2 = new RealParameter(alphaInput);	
		
		// root value
		// Double[] rootValueInput = new Double[] { 2.228585e-40 };
		Double[] rootValueInput = new Double[] { 3.182460e-10 };
		RealParameter rootValue1 = new RealParameter(rootValueInput);
		rootValueInput = new Double[] { 1.021864e-13 };
		RealParameter rootValue2 = new RealParameter(rootValueInput);
		
		// likelihood 1 (condition on rootValue, theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk1 = new OUMVNLikelihoodOneTrait();
		OULk1.initByName("tree", myTree, "sigmasq", sigmasq1, "alpha", alpha1, "optimumManager", optima1,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue1, "eqDist", false);
		lnLk1 = OULk1.calculateLogP();
		
		// likelihood 2 (condition on rootValue, theta_0 = first theta)
		OUMVNLikelihoodOneTrait OULk2 = new OUMVNLikelihoodOneTrait();
				OULk2.initByName("tree", myTree, "sigmasq", sigmasq1, "alpha", alpha1, "optimumManager", optima2,
						"useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue1, "eqDist", false);
		lnLk2 = OULk2.calculateLogP();
		
		// likelihood 3 (rootValue as r.v., theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk3 = new OUMVNLikelihoodOneTrait();
				OULk3.initByName("tree", myTree, "sigmasq", sigmasq2, "alpha", alpha2, "optimumManager", optima3,
						"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue2, "eqDist", true);
		lnLk3 = OULk3.calculateLogP();
		
		// likelihood 4 (rootValue as r.v., theta_0 = first theta)
		OUMVNLikelihoodOneTrait OULk4 = new OUMVNLikelihoodOneTrait();
				OULk4.initByName("tree", myTree, "sigmasq", sigmasq2, "alpha", alpha2, "optimumManager", optima4,
						"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue2, "eqDist", true);
		lnLk4 = OULk4.calculateLogP();
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(9.916106, lnLk1, EPSILON); 
		Assert.assertEquals(9.916106, lnLk2, EPSILON);
		Assert.assertEquals(9.916107, lnLk3, EPSILON);
		Assert.assertEquals(9.916107, lnLk4, EPSILON); 
	}
}
