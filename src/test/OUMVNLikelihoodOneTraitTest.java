package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.ColorManager;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

public class OUMVNLikelihoodOneTraitTest {

	double lnLk1, lnLk2, lnLk3, lnLk4;
	final static double EPSILON = 1e-6;
	
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((sp1[&Regime=1]:1.0,sp2[&Regime=1]:1.0)[&Regime=1]:1.0,sp3[&Regime=2]:2.0)[&Regime=0]:1.0,sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data		
		RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9, 0.0 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
		
		RealParameter colorValues1 = new RealParameter(new Double[] { -4.047373e-16, 4.3, 5.9 }); // thetas
		RealParameter colorValues2 = new RealParameter(new Double[] { 8.128044e-27, 4.3, 5.9 });
		RealParameter colorValues3 = new RealParameter(new Double[] { -1.903330e-16, 4.3, 5.9 });
		RealParameter colorValues4 = new RealParameter(new Double[] { 2.297268e-37, 4.3, 5.9 }); 
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 1, 1, 2, 0, 1, 0, 0 });
							
		ColorManager optima1 = new ColorManager();
		optima1.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues1, "colorAssignments", colorAssignments, "coalCorrection", false);
		ColorManager optima2 = new ColorManager();
		optima2.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues2, "colorAssignments", colorAssignments, "coalCorrection", false);
		ColorManager optima3 = new ColorManager();
		optima3.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues3, "colorAssignments", colorAssignments, "coalCorrection", false);
		ColorManager optima4 = new ColorManager();
		optima4.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues4, "colorAssignments", colorAssignments, "coalCorrection", false);
		
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
		
		// likelihood 1 (condition on rootValue, theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk1 = new OUMVNLikelihoodOneTrait();
		OULk1.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha1, "optimumManager", optima1,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);
		lnLk1 = OULk1.calculateLogP();
		
		// likelihood 2 (condition on rootValue, theta_0 = first theta)
		OUMVNLikelihoodOneTrait OULk2 = new OUMVNLikelihoodOneTrait();
				OULk2.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha1, "optimumManager", optima2,
						"useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);
		lnLk2 = OULk2.calculateLogP();
		
		// likelihood 3 (rootValue as r.v., theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk3 = new OUMVNLikelihoodOneTrait();
				OULk3.initByName("tree", myTree, "sigmasq", sigmasq2, "alpha", alpha2, "optimumManager", optima3,
						"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		lnLk3 = OULk3.calculateLogP();
		
		// likelihood 4 (rootValue as r.v., theta_0 = first theta)
		OUMVNLikelihoodOneTrait OULk4 = new OUMVNLikelihoodOneTrait();
				OULk4.initByName("tree", myTree, "sigmasq", sigmasq2, "alpha", alpha2, "optimumManager", optima4,
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
