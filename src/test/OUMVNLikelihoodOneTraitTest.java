package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.GeneralUtils;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

public class OUMVNLikelihoodOneTraitTest {

	double lnLk1, lnLk2, lnLk3, lnLk4;
	final static double EPSILON = 1e-6;
	
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data		
		String oneTraitValues = "sp1=4.1,sp2= 4.5,sp3=5.9,sp4 =0.0 ";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 1.248328 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 1.734117 };
		RealParameter sigmasq2 = new RealParameter(sigmasqInput);
		
		// alpha
		Double[] alphaInput = new Double[] { 31.20814 };
		RealParameter alpha1 = new RealParameter(alphaInput);
		alphaInput = new Double[] { 43.35287 };
		RealParameter alpha2 = new RealParameter(alphaInput);	
		
		// root value
		Double[] rootValueInput = new Double[] { 2.228585e-40 };
		RealParameter rootValue = new RealParameter(rootValueInput);
		
		// theta
		Double[] thetaInput = { -4.047373e-16, 4.3, 5.9 };
		RealParameter theta1 = new RealParameter(thetaInput);
		thetaInput = new Double[] { 8.128044e-27, 4.3, 5.9 };
		RealParameter theta2 = new RealParameter(thetaInput);
		thetaInput = new Double[] { -1.903330e-16, 4.3, 5.9 };
		RealParameter theta3 = new RealParameter(thetaInput);
		thetaInput = new Double[] { 2.297268e-37, 4.3, 5.9 };
		RealParameter theta4 = new RealParameter(thetaInput);
		
		// likelihood 1 (condition on rootValue, theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk1 = new OUMVNLikelihoodOneTrait();
		OULk1.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha1, "theta", theta1, "nOptima", 3,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);
		lnLk1 = OULk1.calculateLogP();
		
		// likelihood 2 (condition on rootValue, theta_0 = first theta)
		OUMVNLikelihoodOneTrait OULk2 = new OUMVNLikelihoodOneTrait();
				OULk2.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha1, "theta", theta2, "nOptima", 3,
						"useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);
		lnLk2 = OULk2.calculateLogP();
		
		// likelihood 3 (rootValue as r.v., theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk3 = new OUMVNLikelihoodOneTrait();
				OULk3.initByName("tree", myTree, "sigmasq", sigmasq2, "alpha", alpha2, "theta", theta3, "nOptima", 3,
						"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		lnLk3 = OULk3.calculateLogP();
		
		// likelihood 4 (rootValue as r.v., theta_0 = first theta)
		OUMVNLikelihoodOneTrait OULk4 = new OUMVNLikelihoodOneTrait();
				OULk4.initByName("tree", myTree, "sigmasq", sigmasq2, "alpha", alpha2, "theta", theta4, "nOptima", 3,
						"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		lnLk4 = OULk4.calculateLogP();
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(2.148292, lnLk1, EPSILON); 
		Assert.assertEquals(2.148292, lnLk2, EPSILON);
		Assert.assertEquals(2.148292, lnLk3, EPSILON);
		Assert.assertEquals(2.148292, lnLk4, EPSILON); 
	}
}
