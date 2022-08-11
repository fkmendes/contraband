package testdrivers;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import contraband.mvnlikelihood.BMMVNLikelihoodOneTrait;
import contraband.coalescent.CoalCorrection;

import java.util.Arrays;
import java.util.List;

/**
 * @author Fabio K. Mendes
 */

/*
 * Matches testBMMVNLkOneTraitSmallTreeDiffPopSizes
 */
public class BMMVNLikelihoodTestDriver1 {

	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);
		
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
		String spNames = "sp1 sp2 sp3 sp4";
		List<Double> oneTraitValues = Arrays.asList(0.07680552, -0.07201447, -0.03776352, 0.29705797);
		RealParameter oneTraitData = new RealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.006319092 };
		RealParameter sigmaSq1 = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 0.005358022 };
		RealParameter sigmaSq2 = new RealParameter(sigmasqInput);
		
		// root value vector
		Double[] rootValueVectorInput = new Double[] { 0.1 };
		RealParameter rootValue1 = new RealParameter(rootValueVectorInput);
		rootValueVectorInput = new Double[] { 0.09658954 };
		RealParameter rootValue2 = new RealParameter(rootValueVectorInput);

		// likelihood
		BMMVNLikelihoodOneTrait BMLk1 = new BMMVNLikelihoodOneTrait();
		BMLk1.initByName("tree", myTree, "sigmaSq", sigmaSq1, "rootValue", rootValue1, "oneTraitData", oneTraitData, "doCoalCorrection", true, "coalCorrector", coal1);
		BMMVNLikelihoodOneTrait BMLk2 = new BMMVNLikelihoodOneTrait();
		BMLk2.initByName("tree", myTree, "sigmaSq", sigmaSq2, "rootValue", rootValue2, "oneTraitData", oneTraitData, "doCoalCorrection", true, "coalCorrector", coal2);

		System.out.println(BMLk1.calculateLogP()); // -2.1902980
		System.out.println(BMLk2.calculateLogP()); // -2.148631
	}
}
