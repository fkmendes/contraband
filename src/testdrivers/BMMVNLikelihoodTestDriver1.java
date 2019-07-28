package testdrivers;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.BMMVNLikelihoodOneTrait;
import contraband.CoalCorrection;
import contraband.OneValueContTraits;

/*
 * Matches BMMVNLikelihoodTest4
 */
public class BMMVNLikelihoodTestDriver1 {
	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// coal correction
		Double[] smallpopSizesInput = new Double[] { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
		Double[] largePopSizesInput = new Double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
		RealParameter smallPopSizes = new RealParameter(smallpopSizesInput);
		RealParameter largePopSizes = new RealParameter(largePopSizesInput);
		
		CoalCorrection coal1 = new CoalCorrection();
		coal1.initByName("tree", myTree, "popSizes", smallPopSizes);
		CoalCorrection coal2 = new CoalCorrection();
		coal2.initByName("tree", myTree, "popSizes", largePopSizes);
		
		// initializing data
		RealParameter oneTraitValues = new RealParameter(new Double[] { 0.07680552, -0.07201447, -0.03776352, 0.29705797 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.006319092 };
		RealParameter sigmasq1 = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 0.005358022 };
		RealParameter sigmasq2 = new RealParameter(sigmasqInput);
		
		// mean vector
		Double[] meanVectorInput = new Double[] { 0.1 };
		RealParameter mean1 = new RealParameter(meanVectorInput);
		meanVectorInput = new Double[] { 0.09658954 };
		RealParameter mean2 = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNLikelihoodOneTrait BMLk1 = new BMMVNLikelihoodOneTrait();
		BMLk1.initByName("tree", myTree, "sigmasq", sigmasq1, "mean", mean1, "oneTraitData", oneTraitData, "coalCorrector", coal1);
		BMMVNLikelihoodOneTrait BMLk2 = new BMMVNLikelihoodOneTrait();
		BMLk2.initByName("tree", myTree, "sigmasq", sigmasq2, "mean", mean2, "oneTraitData", oneTraitData, "coalCorrector", coal2);
//		BMMVNLikelihoodOneTrait BMLk3 = new BMMVNLikelihoodOneTrait();
//		BMLk3.initByName("tree", myTree, "sigmasq", sigmasq3, "mean", mean, "oneTraitData", oneTraitData);
			
		System.out.println(BMLk1.calculateLogP()); // -2.1902980
		System.out.println(BMLk2.calculateLogP()); // -2.148631
	}
}
