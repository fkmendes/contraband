package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.mvnlikelihood.BMMVNShiftLikelihoodOneTrait;
// import contraband.clock.ColorManager;
import contraband.clock.RateCategoryClockModel;
import contraband.valuewrappers.OneValueContTraits;
import contraband.clock.TreeToVCVMat;

/*
 * This TestDriver matches BMMVNShiftLikelihoodOneTest5. large tree, shift-rate (3 rates) BM
 */
public class BMMVNShiftLikelihoodTestDriver3 {

	public static void main(String[] args) {
		// tree
		String treeStr = "(((((t35:2.336518061,t32:2.336518061):28.95257479,t10:31.28909285):8.654086516,t18:39.94317937):52.28906298,(((t47:31.00652286,t9:31.00652286):50.20634817,(t43:15.06939472,t38:15.06939472):66.14347631):10.61662549,(((((t20:20.94406932,t14:20.94406932):28.09437292,t19:49.03844224):31.16991698,(t24:54.88723469,(((t50:2.534803909,t8:2.534803909):35.18774941,t25:37.72255332):15.41876911,(t12:42.01655137,t5:42.01655137):11.12477106):1.745912255):25.32112453):2.610368667,t37:82.81872788):6.617999642,(t42:81.65977864,(t13:5.88018515,t41:5.88018515):75.77959349):7.776948892):2.392768999):0.4027458181):7.767757656,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0.0;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.006827732, 0.08324919, 0.02754859 });
		
		// first 5 rows are tips
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {
				0, 0, 0, 2, 2, 2, 0, 0, 0, 0,
				2, 1, 2, 1, 1, 1, 2, 2, 1, 0,
				1, 1, 0, 0, 0, 0, 1, 0, 0, 0,
				2, 0, 1, 0, 1, 2, 2, 0, 1, 0,
				1, 0, 1, 0, 2, 2, 0, 0, 2, 0,
				
				0, 0, 0, 0, 0, 0, 2, 2, 2, 2,
				2, 2, 2, 2, 2, 2, 2, 2, 0, 0,
				0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
				1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0
				});
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);
		
		TreeToVCVMat colors = new TreeToVCVMat();
		colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);
		// ColorManager colors = new ColorManager();
		// colors.initByName("tree", myTree, "nTraits", 1, "nColors", 3, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		// initializing data
		RealParameter oneTraitValues = new RealParameter(new Double[] { -1.35423895939605, 0.33827577951559, -1.93240920860189, 3.02404085267233, 0.0465746524280288, -1.73976183419385, 0.484378951244012, 0.65517622307796, 0.75261492151074, -0.24100390220557, 0.947382859472838, 1.98034063758403, 0.933626899580614, -2.03408692528889, -1.09922156803622, -0.0680504653426602, -0.0740826111309686, 1.43932043471363, 1.77514951830302, 0.919758489504323, 1.09574482869026, 0.782136300731067, 0.0962378659531526, -5.33631896191408, 1.4905959236315, -0.0744248296375519, -0.177479586475342, -4.26961572362966, -2.86137408169925, -0.974994208280936, 1.41858116564816, -0.187048179444733, 0.282668872536788, -0.244597214204735, -1.19517123879371, -1.10037443643529, -0.261483237873604, -0.207560694831451, 0.348278374786552, -0.168277479250151, -0.406155392866917, -0.3594972601634, -0.0234868583801822, 1.00637046502527, -0.000469931947963076, -0.329563226090097, 0.043226745691999, 0.739816593725853, -0.238002140615943, -0.214951085346446 });
		String spNames = "t37,t42,t24,t19,t12,t5,t18,t25,t10,t47,t9,t20,t14,t43,t38,t13,t41,t50,t8,t35,t32,t34,t6,t48,t39,t17,t15,t40,t46,t29,t22,t16,t28,t31,t45,t23,t7,t4,t36,t11,t49,t3,t27,t26,t33,t21,t2,t44,t30,t1";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// root value vector
		Double[] rootValueVectorInput = new Double[] { -0.03060231 };
		RealParameter rootValue = new RealParameter(rootValueVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
		System.out.println(BMLk.calculateLogP()); // -65.26436
	}	
}
