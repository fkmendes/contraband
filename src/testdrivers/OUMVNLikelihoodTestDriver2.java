package testdrivers;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

public class OUMVNLikelihoodTestDriver2 {

	private static double lnLk;
	
	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1[&Regime=0]:2.0, sp2[&Regime=0]:1.0)[&Regime=0]:1.0, sp3[&Regime=0]:4.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
				
		// initializing data		
		String oneTraitValues = "sp1=4.1,sp2= 4.5,sp3=5.9,sp4 =0.0 ";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues);
				
		// sigmasq
		Double[] sigmasqInput = new Double[] { 4.003551 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
				
		// alpha
		Double[] alphaInput = new Double[] { 0.7465763 };
		RealParameter alpha = new RealParameter(alphaInput);	
		
		// root value
		Double[] rootValueInput = new Double[] { -33.591241 };
		RealParameter rootValue = new RealParameter(rootValueInput);
				
		// theta
		Double[] thetaInput = { 6.449917 };
		RealParameter theta = new RealParameter(thetaInput);
				
		// likelihood
		OUMVNLikelihoodOneTrait OULk = new OUMVNLikelihoodOneTrait();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "theta", theta, "nOptima", 1,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		
		lnLk = OULk.calculateLogP();
		System.out.println(lnLk);
	}
}
