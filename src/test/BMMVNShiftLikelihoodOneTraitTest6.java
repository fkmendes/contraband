package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.BMMVNShiftLikelihoodOneTrait;
// import contraband.ColorManager;
import contraband.RateCategoryClockModel;
import contraband.OneValueContTraits;
import contraband.TreeToVCVMat;

/*
 * One rate, one trait, on tree with sampled ancestors and fossils (10 extant tips)
 */
public class BMMVNShiftLikelihoodOneTraitTest6 {

	double lnLk;
	final static double EPSILON = 1e-4;
	
	/*
	 * Large tree, with fossils, one-rate BM (on shift-BM class)
	 */
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "((((t16_1:53.23751153,t5_1:53.23751153):8.655770254,((t2_3:8.807060933,t2_2:0):35.79300727,t2_1:0):17.29321358):42.65570341,((((t13_4:37.22097358,t13_3:0):16.65221927,t13_2:0):5.742661592,t13_1:0):32.39904669,((t7_1:29.60227551,t8_1:29.60227551):52.89177712,t28_1:0):9.520848496):12.53408407):67.09599568,((((t15_1:14.38089478,t3_1:0.5275386099):115.9798433,((t1_1:29.18219406,(((t10_3:2.464982927,t10_2:0):1.709987849,t10_1:0):22.8996627,t14_1:2.142260068):42.13119303):28.80221822,t11_1:8.156125409):32.35269336):29.28086375,t19_1:0):7.958804008,((((t4_1:110.1538492,((t6_1:10.73085343,t12_1:16.35390184):26.57130888,t9_1:29.24728433):67.22863851):17.9062182,t20_3:0):5.616707849,t20_2:0):1.720902578,t20_1:0):32.20272799):4.044575033):20.90108373;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.06715078 });
		
		// first 5 rows are tips
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel lsc = new RateCategoryClockModel();
		lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);
		
		TreeToVCVMat colors = new TreeToVCVMat();
		colors.initByName("branchRateModel", lsc, "tree", myTree, "coalCorrection", false);
		// ColorManager colors = new ColorManager();
		// colors.initByName("tree", myTree, "nTraits", 1, "nColors", 1, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
				
		// initializing data	
		RealParameter oneTraitValues = new RealParameter(new Double[] { -5.05744071356103, -1.11052196957346, -0.575427737395457, -0.202953693279079, 1.14969399531796, -1.64202406456342, 1.79417473949251, 0.073438714146757, 0.726726641446067, 3.93444954633643, -0.491682932895839, -1.74464882096195, -0.0780925842656348, 4.81261164374873, 4.49854934867917, 3.04106820281915, -0.803750712288891, 0.523401664662935, -0.456279043339767, 1.2741332640877, 2.87956906600202, 3.91356152320114, 3.59597823070803, -0.389239530235073, -0.63498836297618, -0.789764345711418, -1.20679458892863, -1.48400428008497 });
		String spNames = "t4_1,t12_1,t2_3,t7_1,t13_4,t15_1,t16_1,t5_1,t8_1,t10_3,t6_1,t9_1,t3_1,t14_1,t1_1,t11_1,t2_1,t2_2,t13_1,t13_2,t13_3,t10_1,t10_2,t19_1,t20_1,t20_2,t20_3,t28_1";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// root value vector
		Double[] rootValueVectorInput = new Double[] { -0.4990207 };
		RealParameter rootValue = new RealParameter(rootValueVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
		lnLk = BMLk.calculateLogP();
		System.out.println(lnLk); // -48.54424480901196
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-48.54424, lnLk, EPSILON); 
	}
}
