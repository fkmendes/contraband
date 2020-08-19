package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.clock.ColorManager;
import contraband.mvnlikelihood.OUMVNLikelihoodOneTrait;
import contraband.valuewrappers.OneValueContTraits;

public class OUMVNLikelihoodJiveRealExampleTestDriver {

	private static double lnLk;
	
	public static void main(String[] args) {
		// tree
		String treeStr = "(((((cybotes:11.04802957,(armouri:9.85746743,shrevei:9.85746743):1.190562144):3.31644494,whitemani:14.36447451):2.783755719,(strahmi:12.93828276,longitibialis:12.93828276):4.209947475):5.382653011,marcanoi:22.53088325):21.01029191,(((((quadriocellifer:9.000421692,sagrei:9.000421698):2.997543893,ophiolepis:11.99796559):3.870686412,mestrei:15.868652):1.66332339,(homolechis:6.276287129,jubar:6.276287127):11.25568826):7.19823243,(allogus:15.13139105,(ahli:8.273184746,rubribarbus:8.273184741):6.858206311):9.598816762):18.81096734);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
				
		RealParameter colorValues = new RealParameter(new Double[] { 1.994290, 1.014571 }); // thetas
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {
				0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1,
				1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
							
		ColorManager optima = new ColorManager();
		optima.initByName("nTraits", 1, "nColors", 2, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// initializing data		
		RealParameter oneTraitValues = new RealParameter(new Double[] { 1.83421311376559, 2.0085147653253, 0.800277308189121, 1.32118468496302, 2.04727270131382, 1.81874020073349, 1.10030665760729, 1.15478541494147, 2.10039994621468, 1.90876684524393, 2.307163544304, 2.03620752467502, 1.88733252605389, 0.511196344722762, 1.23993075236991, 0.974319122659916 });
		String spNames = "ahli,allogus,armouri,cybotes,homolechis,jubar,longitibialis,marcanoi,mestrei,ophiolepis,quadriocellifer,rubribarbus,sagrei,shrevei,strahmi,whitemani";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
				
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.109336 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
				
		// alpha
		Double[] alphaInput = new Double[] { 1.321656 };
		RealParameter alpha = new RealParameter(alphaInput);	
		
		// root value
		Double[] rootValueInput = new Double[] { 3.063979e-25 };
		RealParameter rootValue = new RealParameter(rootValueInput);
				
		// likelihood
		OUMVNLikelihoodOneTrait OULk = new OUMVNLikelihoodOneTrait();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "optimumManager", optima,
				"useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);
		
		lnLk = OULk.calculateLogP(); 
		System.out.println(lnLk); // 2.7794249118993495
	}
}
