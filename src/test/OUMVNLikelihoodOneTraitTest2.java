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

public class OUMVNLikelihoodOneTraitTest2 {

	double lnLk1, lnLk2, lnLk3, lnLk4;
	final static double EPSILON = 1e-4;
	
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((sp1:2.0,sp2:1.0):1.0,sp3:4.0):1.0,sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data		
		RealParameter oneTraitValues = new RealParameter(new Double[] { 0.2376494, 0.2950188, 0.8812251, 0.2062229 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
		
		RealParameter colorValues1 = new RealParameter(new Double[] { 1.035041 }); // thetas
		RealParameter colorValues2 = new RealParameter(new Double[] { 0.4152632 });
		RealParameter colorValues3 = new RealParameter(new Double[] { 1.142785e+10 });
		RealParameter colorValues4 = new RealParameter(new Double[] { 0.3497826 }); 
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0 });
		
		RateCategoryClockModel rcc1 = new RateCategoryClockModel();
		rcc1.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues1, "tree", myTree);
		RateCategoryClockModel rcc2 = new RateCategoryClockModel();
		rcc2.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues2, "tree", myTree);
		RateCategoryClockModel rcc3 = new RateCategoryClockModel();
		rcc3.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues3, "tree", myTree);
		RateCategoryClockModel rcc4 = new RateCategoryClockModel();
		rcc4.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues4, "tree", myTree);
		
		TreeToVCVMat optima1 = new TreeToVCVMat();
		optima1.initByName("branchRateModel", rcc1, "tree", myTree, "coalCorrection", false);
		TreeToVCVMat optima2 = new TreeToVCVMat();
		optima2.initByName("branchRateModel", rcc2, "tree", myTree, "coalCorrection", false);
		TreeToVCVMat optima3 = new TreeToVCVMat();
		optima3.initByName("branchRateModel", rcc3, "tree", myTree, "coalCorrection", false);
		TreeToVCVMat optima4 = new TreeToVCVMat();
		optima4.initByName("branchRateModel", rcc4, "tree", myTree, "coalCorrection", false);
		
		/*
		 * Old parameterization using W matrix and ColorManager
		 */
//		ColorManager optima1 = new ColorManager();
//		optima1.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues1, "colorAssignments", colorAssignments, "coalCorrection", false);
//		ColorManager optima2 = new ColorManager();
//		optima2.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues2, "colorAssignments", colorAssignments, "coalCorrection", false);
//		ColorManager optima3 = new ColorManager();
//		optima3.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues3, "colorAssignments", colorAssignments, "coalCorrection", false);
//		ColorManager optima4 = new ColorManager();
//		optima4.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues4, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.02879764 };
		RealParameter sigmasq1 = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 0.09143114 };
		RealParameter sigmasq2 = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 0.009103832 };
		RealParameter sigmasq3 = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 0.01940518 };
		RealParameter sigmasq4 = new RealParameter(sigmasqInput);
		
		// alpha
		Double[] alphaInput = new Double[] { 0.4316411 };
		RealParameter alpha1 = new RealParameter(alphaInput);
		alphaInput = new Double[] { 0.599265 };
		RealParameter alpha2 = new RealParameter(alphaInput);
		alphaInput = new Double[] { 1.818089e-11 };
		RealParameter alpha3 = new RealParameter(alphaInput);
		alphaInput = new Double[] { 2.707329e-11 };
		RealParameter alpha4 = new RealParameter(alphaInput);
		
		// root value
		Double[] rootValueInput = new Double[] { -1.924925 };
		RealParameter rootValue1 = new RealParameter(rootValueInput);
		rootValueInput = new Double[] { -3.726855e-01 };
		RealParameter rootValue2 = new RealParameter(rootValueInput);
		
		// likelihood 1 (condition on rootValue, theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk1 = new OUMVNLikelihoodOneTrait();
		OULk1.initByName("tree", myTree, "sigmasq", sigmasq1, "alpha", alpha1, "optimumManager", optima1,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue1, "eqDist", true);
		lnLk1 = OULk1.calculateLogP();
		
		// likelihood 2 (condition on rootValue, theta_0 = first theta)
		OUMVNLikelihoodOneTrait OULk2 = new OUMVNLikelihoodOneTrait();
				OULk2.initByName("tree", myTree, "sigmasq", sigmasq2, "alpha", alpha2, "optimumManager", optima2,
						"useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue1, "eqDist", true);
		lnLk2 = OULk2.calculateLogP();
		
		// likelihood 3 (rootValue as r.v., theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk3 = new OUMVNLikelihoodOneTrait();
				OULk3.initByName("tree", myTree, "sigmasq", sigmasq3, "alpha", alpha3, "optimumManager", optima3,
						"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue2, "eqDist", false);
		lnLk3 = OULk3.calculateLogP();
		
		// likelihood 4 (rootValue as r.v., theta_0 = first theta)
		OUMVNLikelihoodOneTrait OULk4 = new OUMVNLikelihoodOneTrait();
				OULk4.initByName("tree", myTree, "sigmasq", sigmasq4, "alpha", alpha4, "optimumManager", optima4,
						"useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue2, "eqDist", false);
		lnLk4 = OULk4.calculateLogP();
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(1.171086, lnLk1, EPSILON); 
		Assert.assertEquals(-0.5143806, lnLk2, EPSILON);
		Assert.assertEquals(1.367972, lnLk3, EPSILON);
		Assert.assertEquals(-0.1459122, lnLk4, EPSILON); 
	}
}
