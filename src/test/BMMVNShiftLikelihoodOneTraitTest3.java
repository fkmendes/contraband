package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.BMMVNShiftLikelihoodOneTrait;
import contraband.ColorManager;
import contraband.OneValueContTraits;

public class BMMVNShiftLikelihoodOneTraitTest3 {

	double lnLk, lnLk2;
	final static double EPSILON = 1e-4;
	
	/*
	 * Large tree, with fossils, one-rate BM (on shift-BM class)
	 */
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((((t35:0.1,t32:0.1):0.1,t10:0.1):0.1,t18:0.1):0.1,(((t47:0.1,t9:0.1):0.1,(t43:0.1,t38:0.1):0.1):0.1,(((((t20:0.1,t14:0.1):0.1,t19:0.1):0.1,(t24:0.1,(((t50:0.1,t8:0.1):0.1,t25:0.1):0.1,(t12:0.1,t5:0.1):0.1):0.1):0.1):0.1,t37:0.1):0.1,(t42:0.1,(t13:0.1,t41:0.1):0.1):0.1):0.1):0.1):0.1,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 10.37545 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
		ColorManager colors = new ColorManager();
		colors.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
				
		// initializing data	
		RealParameter oneTraitValues = new RealParameter(new Double[] { -1.98102089083783,0.368221479864884,-2.88206789097561,4.1584950949635,-0.00844451532825974,-2.53652645094459,1.49253360790008,0.860345193251339,2.34223152425472,-0.799642041763574,2.95604288792148,2.77398773288718,1.29414159577677,-6.45163527858711,-3.4910947199795,-0.134187705609486,-0.142843938373703,1.96381042760098,2.43882757691384,2.88182895484916,3.43830441895749,2.47333215100859,0.304330853565934,-6.07220493234317,0.257366287252572,-0.646689981206873,-0.969408565879601,-5.64312827761184,-4.45283272974081,-2.57622923519127,-0.154025403982667,-1.05882534303589,-1.2255223679612,-1.83844113541221,-4.80014501272491,-2.65348584683595,-0.889683582858697,-0.7191201409195,1.03859041594185,-1.54742048206181,-2.30228393056959,-1.58539047472242,-0.522727237241902,-0.531092662298977,-1.53710013936109,-1.86653599613174,-1.50415600947128,-0.806534414120411,-1.20121412215974,-1.12832020868563 });
		String spNames = "t37,t42,t24,t19,t12,t5,t18,t25,t10,t47,t9,t20,t14,t43,t38,t13,t41,t50,t8,t35,t32,t34,t6,t48,t39,t17,t15,t40,t46,t29,t22,t16,t28,t31,t45,t23,t7,t4,t36,t11,t49,t3,t27,t26,t33,t21,t2,t44,t30,t1";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// mean vector
		Double[] meanVectorInput = new Double[] { 0.4877995 };
		RealParameter mean = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "mean", mean, "oneTraitData", oneTraitData);
		lnLk = BMLk.calculateLogP();
		System.out.println(lnLk);
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-168.0481, lnLk, EPSILON); 
	}
}
