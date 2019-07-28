package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.BMMVNLikelihoodOneTrait;
import contraband.CoalCorrection;
import contraband.OneValueContTraits;

/*
 * Ultrametric trees, with different population sizes
 */
public class BMMVNLikelihoodOneTraitTest4 {

	double lnLk1, lnLk2;
	final static double EPSILON = 1e-5;
	
	@Before
	public void setUp() throws Exception {
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
		
		lnLk1 = BMLk1.calculateLogP();
		lnLk2 = BMLk2.calculateLogP();
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(2.190298, lnLk1, EPSILON);
		Assert.assertEquals(2.148630, lnLk2, EPSILON); 
	}
}
