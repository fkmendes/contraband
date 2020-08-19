package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.mvnlikelihood.OUMVNLikelihoodOneTrait;
import contraband.valuewrappers.OneValueContTraits;
import contraband.clock.RateCategoryClockModel;
import contraband.clock.TreeToVCVMat;

public class OUMVNLikelihoodTestDriver3 {

	private static double lnLk;
	
	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:2.0,sp2:1.0):1.0,sp3:4.0):1.0,sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
			
		RealParameter colorValues = new RealParameter(new Double[] { 1.035041 }); // thetas
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel rcc1 = new RateCategoryClockModel();
		rcc1.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);
		
		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc1, "tree", myTree, "coalCorrection", false, "nColors", 1, "colorValues", colorValues, "colorAssignments", colorAssignments);
//		ColorManager optima = new ColorManager();
//		optima.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// initializing data
		RealParameter oneTraitValues = new RealParameter(new Double[] { 0.2376494, 0.2950188, 0.8812251, 0.2062229 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
				
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.02879764 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
				
		// alpha
		Double[] alphaInput = new Double[] { 0.4316411 };
		RealParameter alpha = new RealParameter(alphaInput);	
		
		// root value
		Double[] rootValueInput = new Double[] { -1.924925 };
		RealParameter rootValue = new RealParameter(rootValueInput);
				
		// likelihood
		OUMVNLikelihoodOneTrait OULk = new OUMVNLikelihoodOneTrait();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "optimumManager", optima,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		
		lnLk = OULk.calculateLogP();
		System.out.println(lnLk); // 1.17108619
	}
}
