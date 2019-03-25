package testdrivers;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

public class OUMVNLikelihoodTestDriver {

	private static double lnLk;
	
	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
				
		// initializing data		
		String oneTraitValues = "sp1=4.1,sp2= 4.5,sp3=5.9,sp4 =0.0 ";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues);
				
		// sigmasq
//		Double[] sigmasqInput = new Double[] { 1.248328 };
		Double[] sigmasqInput = new Double[] { 1.734117 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
				
		// alpha
//		Double[] alphaInput = new Double[] { 31.20814 };
		Double[] alphaInput = new Double[] { 43.35287 };
		RealParameter alpha = new RealParameter(alphaInput);	
		
		// root value
//		Double[] rootValueInput = new Double[] { 2.228585e-40 };
//		Double[] rootValueInput = new Double[] { 8.128044e-27 };
		Double[] rootValueInput = new Double[] { 3.348633e-56 };
		RealParameter rootValue = new RealParameter(rootValueInput);
				
		// theta
		Double[] thetaInput = { -4.047373e-16, 4.3, 5.9 };
		RealParameter theta = new RealParameter(thetaInput);
				
		// likelihood
		OUMVNLikelihoodOneTrait OULk = new OUMVNLikelihoodOneTrait();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "theta", theta, "nOptima", 3,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		
		lnLk = OULk.calculateLogP();
		System.out.println(lnLk);
	}
}
