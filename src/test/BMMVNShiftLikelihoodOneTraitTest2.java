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

public class BMMVNShiftLikelihoodOneTraitTest2 {

	double lnLk, lnLk2;
	final static double EPSILON = 1e-5;
	
	/*
	 * Small tree, simple BM. Second likelihood adds root edge. Should match/reduce to BMMVNLikelihoodOneTraitTest
	 */
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((((t35:2.336518061,t32:2.336518061):28.95257479,t10:31.28909285):8.654086516,t18:39.94317937):52.28906298,(((t47:31.00652286,t9:31.00652286):50.20634817,(t43:15.06939472,t38:15.06939472):66.14347631):10.61662549,(((((t20:20.94406932,t14:20.94406932):28.09437292,t19:49.03844224):31.16991698,(t24:54.88723469,(((t50:2.534803909,t8:2.534803909):35.18774941,t25:37.72255332):15.41876911,(t12:42.01655137,t5:42.01655137):11.12477106):1.745912255):25.32112453):2.610368667,t37:82.81872788):6.617999642,(t42:81.65977864,(t13:5.88018515,t41:5.88018515):75.77959349):7.776948892):2.392768999):0.4027458181):7.767757656,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.02813127 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
		ColorManager colors = new ColorManager();
		colors.initByName("nTraits", 1, "maxNColors", 1, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
				
		// initializing data	
		RealParameter oneTraitValues = new RealParameter(new Double[] { 1.44328035263938,3.50655783665801,2.29555775828889,5.27046225284362,0.530549210673224,-0.368179938156466,0.815471685036162,2.36440615897368,0.175573596016938,0.491308438121729,-0.885496203211391,5.25090092164527,4.5030514546132,2.19144454232363,2.71119474862089,3.80879445057103,3.46409119464282,0.847786178303177,1.60806408491409,0.547111053325457,0.463967494705812,3.02072340616352,2.80516460780011,3.79540161037042,4.27779475259434,2.03094267471497,3.13807747280727,0.944275002969514,1.54863412006962,1.46096257283766,0.501600894450101,4.35582751232989,4.53286880080503,2.22222895258926,4.21286791328812,5.90694686097006,1.08900067374343,2.02041235594081,2.35582048390211,4.36319759634288,3.74261911168354,1.26352894420471,2.06627716143412,3.17509471532883,3.28578252973877,2.24446251234181,3.37402888328658,3.67966179487162,4.45836708983263,4.48179955076661 });
		String spNames = "t37,t42,t24,t19,t12,t5,t18,t25,t10,t47,t9,t20,t14,t43,t38,t13,t41,t50,t8,t35,t32,t34,t6,t48,t39,t17,t15,t40,t46,t29,t22,t16,t28,t31,t45,t23,t7,t4,t36,t11,t49,t3,t27,t26,t33,t21,t2,t44,t30,t1";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// mean vector
		Double[] meanVectorInput = new Double[] { 2.405835 };
		RealParameter mean = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "mean", mean, "oneTraitData", oneTraitData);
		lnLk = BMLk.calculateLogP();
		System.out.println(lnLk);
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-75.79013, lnLk, EPSILON); 
	}
}
