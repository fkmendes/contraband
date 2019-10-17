package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.ColorManager;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

public class OUMVNLikelihoodAnolisTestDriver {

	private static double lnLk;
	
	public static void main(String[] args) {
		// tree
		String treeStr = "(((((cybotes:11.04802957,(armouri:9.85746743,shrevei:9.85746743):1.190562144):3.31644494,whitemani:14.36447451):2.783755719,(strahmi:12.93828276,longitibialis:12.93828276):4.209947475):5.382653011,marcanoi:22.53088325):21.01029191,(((((quadriocellifer:9.000421692,sagrei:9.000421698):2.997543893,ophiolepis:11.99796559):3.870686412,mestrei:15.868652):1.66332339,(homolechis:6.276287129,jubar:6.276287127):11.25568826):7.19823243,(allogus:15.13139105,(ahli:8.273184746,rubribarbus:8.273184741):6.858206311):9.598816762):18.81096734);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
			
		RealParameter colorValues = new RealParameter(new Double[] { 57.80504 }); // thetas
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
							
		ColorManager optima = new ColorManager();
		optima.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// initializing data
		RealParameter oneTraitValues = new RealParameter(new Double[] { 56.9933235594816, 57.908923380851, 56.756906958588, 59.5044049185988, 58.0737844282444, 56.7740415733647, 58.252272525633, 58.5358434764243, 58.3521313523052, 57.3562028140739, 59.4100307424803, 58.141977658688, 56.9992157326968, 55.1982328708339, 58.9727991997542, 57.6505778193025 });
		String spNames = "ahli,allogus,armouri,cybotes,homolechis,jubar,longitibialis,marcanoi,mestrei,ophiolepis,quadriocellifer,rubribarbus,sagrei,shrevei,strahmi,whitemani";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
				
		// sigmasq
		Double[] sigmasqInput = new Double[] { 3.255497 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
				
		// alpha
		Double[] alphaInput = new Double[] { 1.404142 };
		RealParameter alpha = new RealParameter(alphaInput);	
		
		// root value
		Double[] rootValueInput = new Double[] { 57.80504 };
		RealParameter rootValue = new RealParameter(rootValueInput);
				
		// likelihood
		OUMVNLikelihoodOneTrait OULk = new OUMVNLikelihoodOneTrait();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "optimumManager", optima,
				"useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);
		
		lnLk = OULk.calculateLogP();
		System.out.println(lnLk);
	}
}
