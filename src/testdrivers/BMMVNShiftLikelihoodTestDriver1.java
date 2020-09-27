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

/**
 * @author Fabio K. Mendes
 */

/*
 * Reduces to testBMMVNLkOneTraitSmallTree
 */
public class BMMVNShiftLikelihoodTestDriver1 {

	public static void main(String[] args) {
		// tree
		String treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);

		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.2704762 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

		TreeToVCVMat colors = new TreeToVCVMat();
		colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);

		// initializing data
		List<Double> oneTraitValues = Arrays.asList(4.1, 4.5, 5.9);
		String spNames = "sp1 sp2 sp3";
		RealParameter oneTraitData = new RealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);
		
		// root value vector
		Double[] rootValueVectorInput = new Double[] { 4.985714 };
		RealParameter rootValue = new RealParameter(rootValueVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
		double lnLk = BMLk.calculateLogP();

		System.out.println(lnLk); // -3.191339
	}	
}
