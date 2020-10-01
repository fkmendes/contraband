package test;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.coalescent.CoalCorrection;
import contraband.mvnlikelihood.BMMVNLikelihoodOneTrait;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

/**
 * @author Fabio K. Mendes
 */

/*
 * This class contains unit tests for BMMVNLikelihoodOneTrait (without and with coalescent correction)
 */
public class BMMVNLikelihoodOneTraitTest {

	final static double EPSILON = 1e-6;

	Tree myTree;
	String treeStr;
	String spNames;
	Double[] sigmasqInput;
	RealParameter sigmaSq;
	Double[] rootValueVectorInput;
	RealParameter rootValue;
	List<Double> oneTraitValues;
	KeyRealParameter oneTraitData;
	BMMVNLikelihoodOneTrait bmLk;

	/*
	 * Small tree, one trait, one rate, without and with root edge (of 1.0).
	 */
	@Test
	public void testBMMVNLkOneTraitSmallTree() {
		// tree
		treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
		myTree = new TreeParser(treeStr, false, false, true, 0);

		/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
		// RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9 });
		// String spNames = "sp1,sp2,sp3";
		// OneValueContTraits oneTraitData = new OneValueContTraits();
		// oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);

		// initializing data
		spNames = "sp1 sp2 sp3";
		oneTraitValues = Arrays.asList(4.1, 4.5, 5.9);
		oneTraitData = new KeyRealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

		// sigmasq
		sigmasqInput = new Double[] { 0.2704762 };
		sigmaSq = new RealParameter(sigmasqInput);

		// root value vector
		rootValueVectorInput = new Double[] { 4.985714 };
		rootValue = new RealParameter(rootValueVectorInput);

		// lnL1
		bmLk = new BMMVNLikelihoodOneTrait();
		bmLk.initByName("tree", myTree, "sigmaSq", sigmaSq, "rootValue", rootValue, "oneTraitData", oneTraitData);
		double lnLk1 = bmLk.calculateLogP(); // no root edge

		bmLk = new BMMVNLikelihoodOneTrait();
		bmLk.initByName("tree", myTree, "sigmaSq", sigmaSq, "rootValue", rootValue, "oneTraitData", oneTraitData, "rootEdgeLength", 1.0);
		double lnLk2 = bmLk.calculateLogP(); // with root edge

		Assert.assertEquals(-3.191339, lnLk1, EPSILON); // no root edge
		Assert.assertEquals(-3.577933, lnLk2, EPSILON); // with root edge
	}

    /*
     * Large ultrametric tree, one trait, one rate, no root edge.
     */
	@Test
	public void testBMMVNLkOneTraitLargeTree() {
		// tree
		treeStr = "((((((t40:5.88018515,t30:5.88018515):31.84236817,(((t43:2.534803909,t25:2.534803909):12.53459081,t1:15.06939472):15.93712814,t4:31.00652286):6.716030459):4.293998055,t7:42.01655137):47.42017615,((t34:39.94317937,t13:39.94317937):42.87554852,t39:82.81872788):6.617999642):2.795514817,(((t19:31.28909285,t50:31.28909285):21.85222958,t9:53.14132243):28.5184562,t26:81.65977864):10.57246371):7.767757656,((t24:69.64568591,((((t5:7.026526759,t27:7.026526759):10.01835474,t45:17.0448815):30.78336706,t8:47.82824856):20.30558317,(((t41:5.397952755,t37:5.397952755):27.52253045,t42:32.9204832):24.003255,t11:56.9237382):11.21009353):1.511854173):12.42349092,(((t16:69.96444119,((t32:24.5178924,t14:24.5178924):37.67637018,(t2:16.45170245,t49:16.45170245):45.74256014):7.7701786):0.3095854321,((t33:30.05535248,((t23:13.30542076,t6:13.30542076):9.449794386,((t48:12.40052663,t21:12.40052663):2.499221124,t44:14.89974775):7.855467399):7.30013733):37.46376453,((t3:46.09553051,t29:46.09553051):5.973201324,t20:52.06873183):15.45038518):2.75490961):6.151101537,((((t28:24.94062793,t46:24.94062793):5.652833062,(t22:27.99970449,(t35:7.370771106,t38:7.370771106):20.62893338):2.593756503):10.27347673,(t17:25.27026498,t10:25.27026498):15.59667274):21.4744838,((t47:6.373725141,t36:6.373725141):45.36629128,(((t15:0.03406798076,t31:0.03406798076):2.516336885,t18:2.550404866):27.30416891,t12:29.85457377):21.88544265):10.60140509):14.08370664):5.644048672):17.93082317):0;";
		myTree = new TreeParser(treeStr, false, false, true, 0);

		/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
		// RealParameter oneTraitValues = new RealParameter(new Double[] { 1.41408472510661, 2.03290942164679, 3.10469256797152, 1.9999136945184, -1.15285873668979, -1.24658888972236, 2.59454370876117, 4.57589687656171, -0.361504996933669, 0.876976899587567, 2.77663728327529, 2.60678865309289, -0.10146053723258, 0.496684400380151, 0.646456479130133, 1.019478423964, 1.02145010112859, -0.228710140607142, 5.07731568044237, 2.05383097646395, -0.0928194844771737, 1.9097451262702, 0.278256939927423, 0.55349407135133, 2.37219286668649, 2.78692166019919, 3.63552330987873, 3.57105666150365, 2.75465573625259, 1.48558479650894, 1.59213516504041, 1.28783779908688, -0.0820595892511218, 1.52897988077741, -0.420502793197038, 1.04372541482933, 1.46053857516283, 0.0322196110491092, 0.28428247250762, 3.38872869386511, 2.59848683743, 2.45004022729315, 1.49883197354575, 0.051896262331336, -0.860501942807566, 1.92657694505103, 2.32348430061041, -0.906837398653968, -1.57544662950472, -1.53788622927357 });
		// String spNames = "t39,t26,t9,t7,t34,t13,t19,t50,t4,t1,t40,t30,t43,t25,t16,t24,t11,t20,t8,t3,t29,t42,t33,t12,t22,t17,t10,t28,t46,t32,t14,t45,t2,t49,t44,t23,t6,t48,t21,t35,t38,t5,t27,t47,t36,t41,t37,t18,t15,t31";
		// OneValueContTraits oneTraitData = new OneValueContTraits();
		// oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);

		// initializing data
		spNames = "t39 t26 t9 t7 t34 t13 t19 t50 t4 t1 t40 t30 t43 t25 t16 t24 t11 t20 t8 t3 t29 t42 t33 t12 t22 t17 t10 t28 t46 t32 t14 t45 t2 t49 t44 t23 t6 t48 t21 t35 t38 t5 t27 t47 t36 t41 t37 t18 t15 t31";
		oneTraitData = new KeyRealParameter();
        oneTraitValues = Arrays.asList(2.6217033, 2.3837527, 0.66385, 2.0554696, 4.6116682, 4.2175155, 1.6890947, 1.1548256, -0.6414763, 2.02788, -0.4061393, 0.2212728, 3.1044734, 2.3983881, 3.4179581, 1.994293, -0.2234107, -0.1731294, -0.9781789, 0.3593945, -1.1369789, 1.7804662, 4.4069651, -0.3154813, 2.3581627, 2.7665761, 2.6999333, 1.051665, 1.5705502, 1.2644854, 2.0151375, 0.0178691, 3.1121187, 2.4489679, 4.461729, 2.8168993, 2.2227894, 6.3508146, 4.6171196, 2.1775336, 2.5118373, -0.2920445, 0.1249836, 1.775187, 1.5983678, 1.4397621, 1.4127476, 3.1335225, 2.7158012, 2.668979);
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

		// sigmasq
		sigmasqInput = new Double[] { 0.03486093 };
		sigmaSq = new RealParameter(sigmasqInput);

		// root value vector
		rootValueVectorInput = new Double[] { 1.859675 };
		rootValue = new RealParameter(rootValueVectorInput);

		BMMVNLikelihoodOneTrait bmLk = new BMMVNLikelihoodOneTrait();
		bmLk.initByName("tree", myTree, "sigmaSq", sigmaSq, "rootValue", rootValue, "oneTraitData", oneTraitData);

		double lnLk = bmLk.calculateLogP(); // with root edge
		Assert.assertEquals(-78.328212, lnLk, EPSILON);
	}

	/*
	 * Large NON-ultrametric tree, one trait, one rate, no root edge, negative root trait value.
	 */
	@Test
	public void testBMMVNLkOneTraitLargeTreeNonUltra() {
		// tree
		String treeStr = "(((((t35:0.1,t32:0.1):0.1,t10:0.1):0.1,t18:0.1):0.1,(((t47:0.1,t9:0.1):0.1,(t43:0.1,t38:0.1):0.1):0.1,(((((t20:0.1,t14:0.1):0.1,t19:0.1):0.1,(t24:0.1,(((t50:0.1,t8:0.1):0.1,t25:0.1):0.1,(t12:0.1,t5:0.1):0.1):0.1):0.1):0.1,t37:0.1):0.1,(t42:0.1,(t13:0.1,t41:0.1):0.1):0.1):0.1):0.1):0.1,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0;";
		myTree = new TreeParser(treeStr, false, false, true, 0);

		/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
		// RealParameter oneTraitValues = new RealParameter(new Double[] { 2.27209315774825,2.35770577479828,2.28619045504381,2.10885751414985,2.21395120993299,2.21616242171277,2.10440705636749,2.23088719352124,2.01612727139994,1.97250353865042,2.04619347343212,1.98626432483887,2.13708092824067,2.0924140883466,2.1149054702513,2.21828980377753,2.21523666408597,2.11038317942733,2.20594065029627,2.17222408711307,2.15525917958223,4.61079655144365,2.05099175361785,2.00896701757015,3.48910077756198,2.57568303611276,2.12434389950229,3.2130625470875,3.07020561053238,1.96918061689753,1.3901244857742,1.45248420753707,1.48390042153334,3.53093468620019,2.88001212857536,1.73916538132083,1.74793606069044,2.15738753957675,1.72286491883736,1.81294836498846,2.03876259543693,1.59727328820493,3.91840500639177,3.23605812584825,3.23116427461224,3.04963684259685,1.99914518841533,3.83711101138963,5.36818679644893,4.67121421317176 });
		// String spNames = "t35,t32,t10,t18,t47,t9,t43,t38,t20,t14,t19,t24,t50,t8,t25,t12,t5,t37,t42,t13,t41,t34,t4,t36,t7,t29,t22,t46,t40,t28,t33,t21,t26,t48,t39,t2,t44,t23,t11,t49,t45,t31,t16,t30,t1,t17,t15,t3,t27,t6";
		// OneValueContTraits oneTraitData = new OneValueContTraits();
		// oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);

		// initializing data
		oneTraitValues = Arrays.asList(2.27209315774825,2.35770577479828,2.28619045504381,2.10885751414985,2.21395120993299,2.21616242171277,2.10440705636749,2.23088719352124,2.01612727139994,1.97250353865042,2.04619347343212,1.98626432483887,2.13708092824067,2.0924140883466,2.1149054702513,2.21828980377753,2.21523666408597,2.11038317942733,2.20594065029627,2.17222408711307,2.15525917958223,4.61079655144365,2.05099175361785,2.00896701757015,3.48910077756198,2.57568303611276,2.12434389950229,3.2130625470875,3.07020561053238,1.96918061689753,1.3901244857742,1.45248420753707,1.48390042153334,3.53093468620019,2.88001212857536,1.73916538132083,1.74793606069044,2.15738753957675,1.72286491883736,1.81294836498846,2.03876259543693,1.59727328820493,3.91840500639177,3.23605812584825,3.23116427461224,3.04963684259685,1.99914518841533,3.83711101138963,5.36818679644893,4.67121421317176);
		spNames = "t35 t32 t10 t18 t47 t9 t43 t38 t20 t14 t19 t24 t50 t8 t25 t12 t5 t37 t42 t13 t41 t34 t4 t36 t7 t29 t22 t46 t40 t28 t33 t21 t26 t48 t39 t2 t44 t23 t11 t49 t45 t31 t16 t30 t1 t17 t15 t3 t27 t6";
		oneTraitData = new KeyRealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

		// sigmasq
		sigmasqInput = new Double[] { 0.01925192 };
		sigmaSq = new RealParameter(sigmasqInput);

		// root value vector
		rootValueVectorInput = new Double[] { 2.182659 };
		rootValue = new RealParameter(rootValueVectorInput);

		BMMVNLikelihoodOneTrait bmLk = new BMMVNLikelihoodOneTrait();
		bmLk.initByName("tree", myTree, "sigmaSq", sigmaSq, "rootValue", rootValue, "oneTraitData", oneTraitData);

		double lnLk = bmLk.calculateLogP(); // with root edge
		Assert.assertEquals(-10.808402, lnLk, EPSILON);
	}

	/*
	 * Small ultrametric tree with different population sizes.
	 */
	@Test
	public void testBMMVNLkOneTraitSmallTreeDiffPopSizes() {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		myTree = new TreeParser(treeStr, false, false, true, 0);

		// coal correction
		Double[] smallpopSizesInput = new Double[] { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
		Double[] largePopSizesInput = new Double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
		RealParameter smallPopSizes = new RealParameter(smallpopSizesInput);
		RealParameter largePopSizes = new RealParameter(largePopSizesInput);

		CoalCorrection coal1 = new CoalCorrection();
		coal1.initByName("tree", myTree, "popSizes", smallPopSizes);
		CoalCorrection coal2 = new CoalCorrection();
		coal2.initByName("tree", myTree, "popSizes", largePopSizes);

		/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
		// RealParameter oneTraitValues = new RealParameter(new Double[] { 0.07680552, -0.07201447, -0.03776352, 0.29705797 });
		// String spNames = "sp1,sp2,sp3,sp4";
		// OneValueContTraits oneTraitData = new OneValueContTraits();
		// oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);

		// initializing data
		spNames = "sp1 sp2 sp3 sp4";
		oneTraitValues = Arrays.asList(0.07680552, -0.07201447, -0.03776352, 0.29705797);
		oneTraitData = new KeyRealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

		// sigmasq
		sigmasqInput = new Double[] { 0.006319092 };
		RealParameter sigmaSq1 = new RealParameter(sigmasqInput);
		sigmasqInput = new Double[] { 0.005358022 };
		RealParameter sigmaSq2 = new RealParameter(sigmasqInput);

		// root value vector
		rootValueVectorInput = new Double[] { 0.1 };
		RealParameter rootValue1 = new RealParameter(rootValueVectorInput);
		rootValueVectorInput = new Double[] { 0.09658954 };
		RealParameter rootValue2 = new RealParameter(rootValueVectorInput);

		// likelihood
		BMMVNLikelihoodOneTrait BMLk1 = new BMMVNLikelihoodOneTrait();
		BMLk1.initByName("tree", myTree, "sigmaSq", sigmaSq1, "rootValue", rootValue1, "oneTraitData", oneTraitData, "doCoalCorrection", true, "coalCorrector", coal1);
		BMMVNLikelihoodOneTrait BMLk2 = new BMMVNLikelihoodOneTrait();
		BMLk2.initByName("tree", myTree, "sigmaSq", sigmaSq2, "rootValue", rootValue2, "oneTraitData", oneTraitData, "doCoalCorrection", true, "coalCorrector", coal2);

		double lnLk1 = BMLk1.calculateLogP();
		double lnLk2 = BMLk2.calculateLogP();

		Assert.assertEquals(2.190298, lnLk1, EPSILON);
		Assert.assertEquals(2.148630, lnLk2, EPSILON);
	}

}
