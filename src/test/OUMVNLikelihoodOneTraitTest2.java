package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.GeneralUtils;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

public class OUMVNLikelihoodOneTraitTest2 {

	double lnLk1, lnLk2, lnLk3, lnLk4;
	final static double EPSILON = 1e-6;
	
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((sp1[&Regime=0]:2.0, sp2[&Regime=0]:1.0)[&Regime=0]:1.0, sp3[&Regime=0]:4.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data		
		RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9, 0.0 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 4.601164 };
		RealParameter sigmasq1 = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 6.841867 };
		RealParameter sigmasq2 = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 4.003551 };
		RealParameter sigmasq3 = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 1.237864 };
		RealParameter sigmasq4 = new RealParameter(sigmasqInput);
		
		// alpha
		Double[] alphaInput = new Double[] { 0.8609833 };
		RealParameter alpha1 = new RealParameter(alphaInput);
		alphaInput = new Double[] { 0.7085376 };
		RealParameter alpha2 = new RealParameter(alphaInput);
		alphaInput = new Double[] { 0.7465763 };
		RealParameter alpha3 = new RealParameter(alphaInput);
		alphaInput = new Double[] { 1.40338e-08 };
		RealParameter alpha4 = new RealParameter(alphaInput);
		
		// root value
		Double[] rootValueInput = new Double[] { -46.965464 };
		RealParameter rootValue1 = new RealParameter(rootValueInput);
		rootValueInput = new Double[] { -33.591241 };
		RealParameter rootValue2 = new RealParameter(rootValueInput);
		
		// theta
		Double[] thetaInput = { 6.201598 };
		RealParameter theta1 = new RealParameter(thetaInput);
		thetaInput = new Double[] { 3.586504 };
		RealParameter theta2 = new RealParameter(thetaInput);
		thetaInput = new Double[] { 6.449917 };
		RealParameter theta3 = new RealParameter(thetaInput);
		thetaInput = new Double[] { 2.792045 };
		RealParameter theta4 = new RealParameter(thetaInput);
		
		// likelihood 1 (condition on rootValue, theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk1 = new OUMVNLikelihoodOneTrait();
		OULk1.initByName("tree", myTree, "sigmasq", sigmasq1, "alpha", alpha1, "theta", theta1, "nOptima", 1,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue1, "eqDist", true);
		lnLk1 = OULk1.calculateLogP();
		
		// likelihood 2 (condition on rootValue, theta_0 = first theta)
		OUMVNLikelihoodOneTrait OULk2 = new OUMVNLikelihoodOneTrait();
				OULk2.initByName("tree", myTree, "sigmasq", sigmasq2, "alpha", alpha2, "theta", theta2, "nOptima", 1,
						"useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue1, "eqDist", true);
		lnLk2 = OULk2.calculateLogP();
		
		// likelihood 3 (rootValue as r.v., theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk3 = new OUMVNLikelihoodOneTrait();
				OULk3.initByName("tree", myTree, "sigmasq", sigmasq3, "alpha", alpha3, "theta", theta3, "nOptima", 1,
						"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue2, "eqDist", false);
		lnLk3 = OULk3.calculateLogP();
		
		// likelihood 4 (rootValue as r.v., theta_0 = first theta)
		OUMVNLikelihoodOneTrait OULk4 = new OUMVNLikelihoodOneTrait();
				OULk4.initByName("tree", myTree, "sigmasq", sigmasq4, "alpha", alpha4, "theta", theta4, "nOptima", 1,
						"useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue2, "eqDist", false);
		lnLk4 = OULk4.calculateLogP();
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-7.63854, lnLk1, EPSILON); 
		Assert.assertEquals(-8.817273, lnLk2, EPSILON);
		Assert.assertEquals(-7.630117, lnLk3, EPSILON);
		Assert.assertEquals(-8.457486, lnLk4, EPSILON); 
	}
}
