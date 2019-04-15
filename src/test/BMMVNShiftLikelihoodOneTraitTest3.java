package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.BMMVNLikelihoodOneTrait;
import contraband.BMMVNShiftLikelihoodOneTrait;
import contraband.ColorManager;
import contraband.OneValueContTraits;

public class BMMVNShiftLikelihoodOneTraitTest3 {

	double lnLk, lnLk2;
	final static double EPSILON = 1e-5;
	
	/*
	 * Small tree, simple BM. Second likelihood adds root edge. Should match/reduce to BMMVNLikelihoodOneTraitTest
	 */
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((((t35:0.1,t32:0.1):0.1,t10:0.1):0.1,t18:0.1):0.1,(((t47:0.1,t9:0.1):0.1,(t43:0.1,t38:0.1):0.1):0.1,(((((t20:0.1,t14:0.1):0.1,t19:0.1):0.1,(t24:0.1,(((t50:0.1,t8:0.1):0.1,t25:0.1):0.1,(t12:0.1,t5:0.1):0.1):0.1):0.1):0.1,t37:0.1):0.1,(t42:0.1,(t13:0.1,t41:0.1):0.1):0.1):0.1):0.1):0.1,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.01925192 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
		ColorManager colors = new ColorManager();
		colors.initByName("tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
				
		// initializing data	
		RealParameter oneTraitValues = new RealParameter(new Double[] { 2.27209315774825,2.35770577479828,2.28619045504381,2.10885751414985,2.21395120993299,2.21616242171277,2.10440705636749,2.23088719352124,2.01612727139994,1.97250353865042,2.04619347343212,1.98626432483887,2.13708092824067,2.0924140883466,2.1149054702513,2.21828980377753,2.21523666408597,2.11038317942733,2.20594065029627,2.17222408711307,2.15525917958223,4.61079655144365,2.05099175361785,2.00896701757015,3.48910077756198,2.57568303611276,2.12434389950229,3.2130625470875,3.07020561053238,1.96918061689753,1.3901244857742,1.45248420753707,1.48390042153334,3.53093468620019,2.88001212857536,1.73916538132083,1.74793606069044,2.15738753957675,1.72286491883736,1.81294836498846,2.03876259543693,1.59727328820493,3.91840500639177,3.23605812584825,3.23116427461224,3.04963684259685,1.99914518841533,3.83711101138963,5.36818679644893,4.67121421317176 });
		String spNames = "t35,t32,t10,t18,t47,t9,t43,t38,t20,t14,t19,t24,t50,t8,t25,t12,t5,t37,t42,t13,t41,t34,t4,t36,t7,t29,t22,t46,t40,t28,t33,t21,t26,t48,t39,t2,t44,t23,t11,t49,t45,t31,t16,t30,t1,t17,t15,t3,t27,t6";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// mean vector
		Double[] meanVectorInput = new Double[] { 2.182659 };
		RealParameter mean = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "mean", mean, "oneTraitData", oneTraitData);
		lnLk = BMLk.calculateLogP();
		System.out.println(lnLk);
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-10.8084, lnLk, EPSILON); 
	}
}
