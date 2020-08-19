package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import contraband.valuewrappers.OneValueContTraits;
import contraband.otherlikelihood.WNLikelihood;

public class WNLikelihoodOneTraitTest {

	double lnLk;
	final static double EPSILON = 1e-6;
	
	@Before
	public void setUp() throws Exception {
		// initializing data
		Double samples[] = new Double[] { 1.5952808, 0.3295078, -0.8204684 };
		RealParameter oneTraitValues = new RealParameter(samples);
		String spNames = "sp1,sp2,sp3";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// double logVar = Math.log(0.9733856);
		// System.out.println(logVar);
		RealParameter colorValues = new RealParameter(new Double[] { 0.9733856, 0.0, 1.0 });	
		RealParameter means = new RealParameter(new Double[] { 0.3681067, 0.1, 0.2 });
		IntegerParameter normalAssignments = new IntegerParameter(new Integer[] { 0, 0, 0 });
		
		WNLikelihood WNLk = new WNLikelihood();
		WNLk.initByName("oneTraitData", oneTraitData, "sigmaSqs", colorValues, "mus", means, "normalAssignments", normalAssignments);
		lnLk = WNLk.calculateLogP();
		System.out.println(lnLk);
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-4.216353, lnLk, EPSILON);
	}
}
