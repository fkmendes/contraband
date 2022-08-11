package testdrivers;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import contraband.mvnlikelihood.OUMVNLikelihoodOneTrait;
import contraband.clock.RateCategoryClockModel;
import contraband.clock.TreeToVCVMat;

import java.util.Arrays;
import java.util.List;

/*
 * Matches testOUMVNLkOneTraitSmallTree3optRandomRVEstimateRV
 */
public class OUMVNLikelihoodTestDriver1 {

	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);

		// thetas
		RealParameter colorValues = new RealParameter(new Double[]{ 0.206222932117995, 0.26633408087427, 0.88122539543514 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[]{ 1, 1, 2, 0, 1, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

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
		Double[] sigmasqInput = new Double[] { 0.008287661 };
		RealParameter sigmaSq = new RealParameter(sigmasqInput);
				
		// alpha
		Double[] alphaInput = new Double[] { 10.07163 };
		RealParameter alpha = new RealParameter(alphaInput);	
		
		// root value
		Double[] rootValueInput = new Double[] { 1.021864e-13 };
		RealParameter rootValue = new RealParameter(rootValueInput);
				
		// likelihood
		OUMVNLikelihoodOneTrait ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree, "sigmaSq", sigmaSq, "alpha", alpha, "optimumManager", optima, "useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true); // rootValue here is irrelevant
		double lnLk = ouLk.calculateLogP();
		
		lnLk = ouLk.calculateLogP();
		System.out.println(lnLk); // 9.916106821947409
	}
}
