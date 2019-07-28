package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;
import contraband.RateCategoryClockModel;
import contraband.TreeToVCVMat;

/*
 * Equivalent to first likelihood test inside OUMVNLikelihoodOneTraitTest2.java
 */
public class OUMVNLikelihoodTestDriver1 {

	private static double lnLk;
	
	public static void main(String[] args) {
		// tree
		// String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
				
		RealParameter colorValues = new RealParameter(new Double[] { -4.047373e-16, 4.3, 5.9 }); // thetas
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 1, 1, 2, 0, 1, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);
		
		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false, "nColors", 3, "colorValues", colorValues, "colorAssignments", colorAssignments);
		// ColorManager optima = new ColorManager();
		// optima.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// initializing data		
		RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9, 0.0 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
				
		// sigmasq
		Double[] sigmasqInput = new Double[] { 1.734117 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
				
		// alpha
		Double[] alphaInput = new Double[] { 43.35287 };
		RealParameter alpha = new RealParameter(alphaInput);	
		
		// root value
		Double[] rootValueInput = new Double[] { 3.348633e-56 };
		RealParameter rootValue = new RealParameter(rootValueInput);
				
		// likelihood
		OUMVNLikelihoodOneTrait OULk = new OUMVNLikelihoodOneTrait();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "optimumManager", optima,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		
		lnLk = OULk.calculateLogP(); 
		System.out.println(lnLk); // 2.1482918780359883
	}
}
