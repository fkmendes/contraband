package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.mvnlikelihood.OUMVNLikelihoodOneTrait;
import contraband.valuewrappers.OneValueContTraits;
import contraband.clock.RateCategoryClockModel;
import contraband.clock.TreeToVCVMat;

/*
 * Equivalent to second likelihood test inside OUMVNLikelihoodOneTraitTest2.java
 */
public class OUMVNLikelihoodTestDriver2 {

	private static double lnLk;
	
	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
				
		RealParameter colorValues = new RealParameter(new Double[] { 0.4152632 }); // thetas
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);
		
		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false, "nColors", 1, "colorValues", colorValues, "colorAssignments", colorAssignments);
		// ColorManager optima = new ColorManager();
		// optima.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// initializing data		
		RealParameter oneTraitValues = new RealParameter(new Double[] { 0.237649365136715, 0.295018750722361, 0.881225138279161, 0.206222932069516 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
				
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.09143114 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
				
		// alpha
		Double[] alphaInput = new Double[] { 0.599265 };
		RealParameter alpha = new RealParameter(alphaInput);	
		
		// root value
		Double[] rootValueInput = new Double[] { 3.348633e-56 };
		RealParameter rootValue = new RealParameter(rootValueInput);
				
		// likelihood
		OUMVNLikelihoodOneTrait OULk = new OUMVNLikelihoodOneTrait();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "optimumManager", optima,
				"useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		
		lnLk = OULk.calculateLogP(); 
		System.out.println(lnLk); // 2.1482918780359883
	}
}
