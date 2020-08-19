package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.mvnlikelihood.BMMVNShiftLikelihoodOneTrait;
// import contraband.clock.ColorManager;
import contraband.clock.RateCategoryClockModel;
import contraband.valuewrappers.OneValueContTraits;
import contraband.clock.TreeToVCVMat;

/*
 * This testdriver runs the BMMVNShift class with a single rate for the whole tree (matches/reduces to BMMVNLikelihoodTestDriver)
 */
public class BMMVNShiftLikelihoodTestDriver1 {

	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:1.0, sp2:1.0):1.0, sp3:2.0):1.0, sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 1.4822794118 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel lsc = new RateCategoryClockModel();
		lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);
		
		TreeToVCVMat colors = new TreeToVCVMat();
		colors.initByName("branchRateModel", lsc, "tree", myTree, "coalCorrection", false);
		// ColorManager colors = new ColorManager();
		// colors.initByName("tree", myTree, "nTraits", 1, "nColors", 1, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// initializing data
		RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9, 0.0 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// root value vector
		Double[] rootValueVectorInput = new Double[] { 3.079142 };
		RealParameter rootValue = new RealParameter(rootValueVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
		System.out.println(BMLk.calculateLogP()); // -8.29469706
	}	
}
