package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.mvnlikelihood.OUMVNLikelihoodOneTrait;
import contraband.clock.RateCategoryClockModel;
import contraband.clock.TreeToVCVMat;

import java.util.Arrays;
import java.util.List;

/*
 * Matches testOUMVNLkOneTraitSmallNonUltraTree3optRandomRVEstimateRV
 */
public class OUMVNLikelihoodTestDriver2 {
	
	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:2.0,sp2:1.0):1.0,sp3:4.0):1.0,sp4:3.0);";
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);

		// thetas
		RealParameter colorValues = new RealParameter(new Double[]{ 0.4152632 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);
		// ColorManager optima = new ColorManager();
		// optima.initByName("nTraits", 1, "nColors", 3, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// initializing data		
		String spNames = "sp1 sp2 sp3 sp4";
		List<Double> oneTraitValues = Arrays.asList(0.237649365136715, 0.295018750722361, 0.881225138279161, 0.206222932069516);
		RealParameter oneTraitData = new RealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);
				
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.09143114 };
		RealParameter sigmaSq = new RealParameter(sigmasqInput);
				
		// alpha
		Double[] alphaInput = new Double[] { 0.599265 };
		RealParameter alpha = new RealParameter(alphaInput);	
		
		// root value
		Double[] rootValueInput = new Double[] { 1.924925-56 };
		RealParameter rootValue = new RealParameter(rootValueInput);
				
		// likelihood
		OUMVNLikelihoodOneTrait ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree, "sigmaSq", sigmaSq, "alpha", alpha, "optimumManager", optima, "useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		double lnLk = ouLk.calculateLogP();

		System.out.println(lnLk); // -0.5143806
	}
}
