package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;
import contraband.RateCategoryClockModel;
import contraband.TreeToVCVMat;

/*
 * Small non-ultrametric tree with sampled ancestors and fossil tips, two optima
 */
public class OUMVNLikelihoodTestDriver4 {

	private static double lnLk;

	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:2.0,sp2:1.0):1.0,sp3:0.0):1.0,sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);

		// initializing data
		RealParameter oneTraitValues = new RealParameter(new Double[] { 0.267347579902117, 0.322440128331226, 0.584766339144924, 0.211885389593768 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);

		RealParameter colorValues = new RealParameter(new Double[] { 0.05779027, 0.19382641 }); // thetas
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 1, 1, 1, 0, 1, 1, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 2, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);
		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);
//		ColorManager optima = new ColorManager();
//		optima.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);

		// sigmasq
		Double[] sigmasqInput = new Double[] { 1.11256e-08 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);

		// alpha
		Double[] alphaInput = new Double[] { 0.5564338 };
		RealParameter alpha = new RealParameter(alphaInput);

		// root value
		Double[] rootValueInput = new Double[] { 0.87579783 };
		RealParameter rootValue = new RealParameter(rootValueInput);

		// likelihood 1 (condition on rootValue, theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk = new OUMVNLikelihoodOneTrait();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "optimumManager", optima,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);

		lnLk = OULk.calculateLogP();
		System.out.println("LnLk=" + lnLk);
	}
}
