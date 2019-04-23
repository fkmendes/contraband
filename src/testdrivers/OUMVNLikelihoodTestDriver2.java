package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.ColorManager;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OUMVNLikelihoodOneTrait_NoTreeSampling;
import contraband.OneValueContTraits;

public class OUMVNLikelihoodTestDriver2 {

	private static double lnLk;
	
	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:2.0,sp2:1.0):1.0,sp3:4.0):1.0,sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
			
		RealParameter colorValues = new RealParameter(new Double[] { 6.449917 }); // thetas
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0 });
							
		ColorManager optima = new ColorManager();
		optima.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// initializing data
		RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9, 0.0 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
				
		// sigmasq
		Double[] sigmasqInput = new Double[] { 4.003551 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
				
		// alpha
		Double[] alphaInput = new Double[] { 0.7465763 };
		RealParameter alpha = new RealParameter(alphaInput);	
		
		// root value
		Double[] rootValueInput = new Double[] { -33.591241 };
		RealParameter rootValue = new RealParameter(rootValueInput);
				
		// likelihood
		OUMVNLikelihoodOneTrait OULk = new OUMVNLikelihoodOneTrait();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "optimumManager", optima,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		
		lnLk = OULk.calculateLogP();
		System.out.println(lnLk);
	}
}
