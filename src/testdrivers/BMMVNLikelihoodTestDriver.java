package testdrivers;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.BMMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

public class BMMVNLikelihoodTestDriver {

	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data
		RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9, 0.0 });
		String spNames = "sp1,sp2,sp3,sp4";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 1.4822794118 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
		
		// mean vector
		Double[] meanVectorInput = new Double[] { 3.079142 };
		RealParameter mean = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNLikelihoodOneTrait BMLk = new BMMVNLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "sigmasq", sigmasq, "mean", mean, "oneTraitData", oneTraitData);
		System.out.println(BMLk.calculateLogP());
	}
}
