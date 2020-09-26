package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.mvnlikelihood.BMMVNShiftLikelihoodOneTrait;
import contraband.clock.RateCategoryClockModel;
import contraband.clock.TreeToVCVMat;

import java.util.Arrays;
import java.util.List;

/*
 * Matches testBMMVNShiftLkOneTraitSmallTreeTwoRates
 */
public class BMMVNShiftLikelihoodTestDriver2 {

	public static void main(String[] args) {
		// tree
		String treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);

		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.05057867, 3.360241 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 1, 1, 1 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 2, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

		TreeToVCVMat colors = new TreeToVCVMat();
		colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);

		// initializing data
		List<Double> oneTraitValues = Arrays.asList(-2.53718502574816, -2.85562629168723, 1.79661600241838);
		String spNames = "sp1 sp2 sp3";
		RealParameter oneTraitData = new RealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);
		
		// root value vector
		Double[] rootValueVectorInput = new Double[] { -1.191236 };
		RealParameter rootValue = new RealParameter(rootValueVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait bmLk = new BMMVNShiftLikelihoodOneTrait();
		bmLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
		double lnLk = bmLk.calculateLogP();

		System.out.println(lnLk); // -4.673609033125993
	}	
}
