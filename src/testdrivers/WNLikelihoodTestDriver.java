package testdrivers;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import contraband.otherlikelihood.WNLikelihoodOneTrait;

import java.util.Arrays;
import java.util.List;

/**
 * @author Fabio K. Mendes
 */

/*
 * This testdriver runs the WNLikelihood class with a single rate for the whole tree
 */
public class WNLikelihoodTestDriver {

	public static void main(String[] args) {

		/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
		// Double samples[] = new Double[] { 1.5952808, 0.3295078, -0.8204684 };
		// RealParameter oneTraitValues = new RealParameter(samples);
		// String spNames = "sp1,sp2,sp3";
		// OneValueContTraits oneTraitData = new OneValueContTraits();
		// oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);

		// initializing data
		String spNames = "sp1 sp2 sp3";
		List<Double> oneTraitValues = Arrays.asList(1.5952808, 0.3295078, -0.8204684);
		RealParameter oneTraitData = new RealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

		RealParameter sigmaSqs = new RealParameter(new Double[] { 0.9733856, 0.0, 1.0 });
		RealParameter mus = new RealParameter(new Double[] { 0.3681067, 0.1, 0.2 });
		IntegerParameter normalAssignments = new IntegerParameter(new Integer[] { 0, 0, 0 });

		WNLikelihoodOneTrait WNLk = new WNLikelihoodOneTrait();
		WNLk.initByName("oneTraitData", oneTraitData, "sigmaSqs", sigmaSqs, "mus", mus, "normalAssignments", normalAssignments);
		double lnLk = WNLk.calculateLogP();

		System.out.println(lnLk); // -4.216353195954754
	}	
}
