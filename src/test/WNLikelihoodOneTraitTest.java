package test;

import org.junit.Assert;
import org.junit.Test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import contraband.otherlikelihood.WNLikelihoodOneTrait;
import outercore.parameter.KeyEnhancedRealParameter;

import java.util.Arrays;
import java.util.List;

/**
 * @author Fabio K. Mendes
 */

public class WNLikelihoodOneTraitTest {

	final static double EPSILON = 1e-6;
	
	@Test
	public void testWNLkOneTrait() {

		/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
		// Double samples[] = new Double[] { 1.5952808, 0.3295078, -0.8204684 };
		// RealParameter oneTraitValues = new RealParameter(samples);
		// String spNames = "sp1,sp2,sp3";
		// OneValueContTraits oneTraitData = new OneValueContTraits();
		// oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);

		// initializing data
		String spNames = "sp1 sp2 sp3";
		List<Double> oneTraitValues = Arrays.asList(1.5952808, 0.3295078, -0.8204684);
		KeyEnhancedRealParameter oneTraitData = new KeyEnhancedRealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

		RealParameter sigmaSqs = new RealParameter(new Double[] { 0.9733856, 0.0, 1.0 });
		RealParameter mus = new RealParameter(new Double[] { 0.3681067, 0.1, 0.2 });
		IntegerParameter normalAssignments = new IntegerParameter(new Integer[] { 0, 0, 0 });
		
		WNLikelihoodOneTrait wnLk = new WNLikelihoodOneTrait();
		wnLk.initByName("oneTraitData", oneTraitData, "sigmaSqs", sigmaSqs, "mus", mus, "normalAssignments", normalAssignments);
		double lnLk = wnLk.calculateLogP();

		Assert.assertEquals(-4.216353, lnLk, EPSILON);
	}
}
