package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.BMMVNShiftLikelihoodOneTrait;
import contraband.ColorManager;
import contraband.RateCategoryClockModel;
import contraband.OneValueContTraits;
import contraband.TreeToVCVMat;

public class BMMVNShiftLikelihoodOneTraitTest4 {

	double lnLk;
	final static double EPSILON = 1e-4;
	
	/*
	 * Small tree with no fossils, two rates BM (on shift-BM class)
	 */
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.05057867, 3.360241 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 1, 1, 1 });
		RateCategoryClockModel lsc = new RateCategoryClockModel();
		lsc.initByName("nCat", 2, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);
		
		TreeToVCVMat colors = new TreeToVCVMat();
		colors.initByName("branchRateModel", lsc, "tree", myTree, "coalCorrection", false);		
		// ColorManager colors = new ColorManager();
		// colors.initByName("tree", myTree, "nTraits", 1, "nColors", 2, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
				
		// initializing data	
		RealParameter oneTraitValues = new RealParameter(new Double[] { -2.53718502574816, -2.85562629168723, 1.79661600241838 });
		String spNames = "sp1,sp2,sp3";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// mean vector
		Double[] meanVectorInput = new Double[] { -1.191236 };
		RealParameter mean = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "mean", mean, "oneTraitData", oneTraitData);
		lnLk = BMLk.calculateLogP();
		System.out.println(lnLk); // -4.673609  
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-4.673609  , lnLk, EPSILON); 
	}
}
