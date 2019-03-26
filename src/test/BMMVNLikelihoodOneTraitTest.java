package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.BMMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

public class BMMVNLikelihoodOneTraitTest {

	double lnLk, lnLk2;
	final static double EPSILON = 1e-6;
	
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "((sp1:1.0, sp2:1.0):1.0,sp3:2.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data		
		String oneTraitValues = "sp1=4.1,sp2= 4.5, sp3=5.9";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.270476 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
		
		// mean vector
		Double[] meanVectorInput = new Double[] { 4.9857143 };
		RealParameter mean = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNLikelihoodOneTrait BMLk = new BMMVNLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "sigmasq", sigmasq, "mean", mean, "oneTraitData", oneTraitData);
		lnLk = BMLk.calculateLogP(); // no root edge
		
		BMMVNLikelihoodOneTrait BMLk2 = new BMMVNLikelihoodOneTrait();
		BMLk2.initByName("tree", myTree, "sigmasq", sigmasq, "mean", mean, "oneTraitData", oneTraitData, "rootEdgeLength", 1.0);
		lnLk2 = BMLk2.calculateLogP(); // with root edge
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-3.191339, lnLk, EPSILON);
		Assert.assertEquals(-3.577933, lnLk2, EPSILON);
	}
}
