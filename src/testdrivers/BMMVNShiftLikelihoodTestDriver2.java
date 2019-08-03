package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.BMMVNShiftLikelihoodOneTrait;
import contraband.RateCategoryClockModel;
import contraband.OneValueContTraits;
import contraband.TreeToVCVMat;

/*
 * This TestDriver matches BMMVNShiftLikelihoodOneTest4. Small tree, shift-rate (2 rates) BM
 */
public class BMMVNShiftLikelihoodTestDriver2 {

	public static void main(String[] args) {
		// tree
		String treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.1160941, 0.01914707 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 1, 1, 0, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 2, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);
		
		TreeToVCVMat colors = new TreeToVCVMat();
		colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);
		// ColorManager colors = new ColorManager();
		// colors.initByName("tree", myTree, "nTraits", 1, "nColors", 2, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// initializing data
		RealParameter oneTraitValues = new RealParameter(new Double[] { -0.9291812, -0.7312343, -1.6712572 });
		String spNames = "sp1,sp2,sp3";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// root value vector
		Double[] rootValueVectorInput = new Double[] { -1.125558 };
		RealParameter rootValue = new RealParameter(rootValueVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
		System.out.println(BMLk.calculateLogP()); // -0.8583676
	}	
}
