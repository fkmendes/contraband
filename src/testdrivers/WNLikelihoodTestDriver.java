package testdrivers;

import beast.core.parameter.RealParameter;
import contraband.MVNUtils;
import contraband.OneValueContTraits;
import contraband.WNLikelihood;

/*
 * This testdriver runs the WNLikelihood class with a single rate for the whole tree
 */
public class WNLikelihoodTestDriver {

	public static void main(String[] args) {	
		Double samples[] = new Double[] { 1.5952808, 0.3295078, -0.8204684 };
		Double sigmaSqs[] = new Double[] { 0.9733856, 0.9733856, 0.9733856 };
		Double mus[] = new Double[] { 0.3681067, 0.3681067, 0.3681067 };
		System.out.println(MVNUtils.getSampleMultipleNormalLogLk(samples, mus, sigmaSqs));
		
		// initializing data
		RealParameter oneTraitValues = new RealParameter(samples);
		String spNames = "sp1,sp2,sp3";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		RealParameter colorValues = new RealParameter(sigmaSqs);
		RealParameter means = new RealParameter(mus);
		
		WNLikelihood WNLk = new WNLikelihood();
		WNLk.initByName("oneTraitData", oneTraitData, "sigmaSqs", colorValues, "mus", means);
		System.out.println(WNLk.calculateLogP());
	}	
}
