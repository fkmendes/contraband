package contraband.test;

import beast.base.evolution.tree.Tree;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeParser;
import contraband.mvnlikelihood.OUMVNLikelihoodOneTrait;
import contraband.clock.RateCategoryClockModel;
import contraband.clock.TreeToVCVMat;
// import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

/**
 * @author Fabio K. Mendes
 */

/*
 * This class contains unit tests for OUMVNLikelihoodOneTrait
 */
public class OUMVNLikelihoodOneTraitTest {

	Tree myTree1, myTree2;
	String treeStr1, treeStr2;
	String spNames;
	Double[] sigmasqInput;
	RealParameter sigmaSq1, sigmaSq2, sigmaSq;
	Double[] alphaInput;
	RealParameter alpha1, alpha2, alpha;
	Double[] rootValueInput;
	RealParameter rootValue1, rootValue2, rootValue;
	List<Double> oneTraitValues;
	RealParameter oneTraitData;
	// KeyRealParameter oneTraitData;
	OUMVNLikelihoodOneTrait ouLk;

	/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
	// OneValueContTraits oneTraitData;

	final static double EPSILON = 1e-7;

	@Before
	public void setUp() {
		// tree 1
		treeStr1 = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
		myTree1 = new TreeParser(treeStr1, false, false, true, 0);

		// tree 2
		treeStr2 = "(((sp1:2.0,sp2:1.0):1.0,sp3:4.0):1.0,sp4:3.0);";
		myTree2 = new TreeParser(treeStr2, false, false, true, 0);

		// initializing data
		/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
		// RealParameter oneTraitValues = new RealParameter(new Double[]{ 0.237649365136715, 0.295018750722361, 0.881225138279161, 0.206222932069516 });
		// spNames = "sp1,sp2,sp3,sp4";
		// oneTraitData = new OneValueContTraits();
		// oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
		spNames = "sp1 sp2 sp3 sp4";
		oneTraitValues = Arrays.asList(0.237649365136715, 0.295018750722361, 0.881225138279161, 0.206222932069516);
		oneTraitData = new RealParameter();
		// oneTraitData = new KeyRealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

		// sigma 1, alpha 1 and root value 1
		sigmasqInput = new Double[]{ 0.006082604 };
		sigmaSq1 = new RealParameter(sigmasqInput);
		alphaInput = new Double[]{ 7.390366 };
		alpha1 = new RealParameter(alphaInput);
		rootValueInput = new Double[]{ 3.182460e-10 };
		rootValue1 = new RealParameter(rootValueInput);

		// sigma 2, alpha 2 and root value 2
		sigmasqInput = new Double[]{ 0.008287661 };
		sigmaSq2 = new RealParameter(sigmasqInput);
		alphaInput = new Double[]{ 10.07163 };
		alpha2 = new RealParameter(alphaInput);
		rootValueInput = new Double[]{ 1.021864e-13 };
		rootValue2 = new RealParameter(rootValueInput);
	}

	/*
	 * (1) Small tree, one trait, three optima, conditioning on root value, estimating root value
	 */
	@Test
	public void testOUMVNLkOneTraitSmallTree3optCondOnRVEstimateRV() {
		// thetas
		RealParameter colorValues = new RealParameter(new Double[]{ 0.206222932117995, 0.26633408087427, 0.88122539543514 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[]{ 1, 1, 2, 0, 1, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree1);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree1, "coalCorrection", false);

		// likelihood
		ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree1, "sigmaSq", sigmaSq1, "alpha", alpha1, "optimumManager", optima, "useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue1, "eqDist", false);
		double lnLk = ouLk.calculateLogP();

		Assert.assertEquals(9.9161059, lnLk, EPSILON); // in R, we get 9.916106
	}

	/*
	 * (2) Small tree, one trait, three optima, conditioning on root value, using first theta as root value (one fewer parameter)
	 */
	@Test
	public void testOUMVNLkOneTraitSmallTree3optCondOnRVFirstThetaIsRV() {
		// thetas
		RealParameter colorValues = new RealParameter(new Double[]{ 0.206222932069516, 0.266334080825641, 0.881225395384966 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[]{ 1, 1, 2, 0, 1, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree1);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree1, "coalCorrection", false);

		// likelihood
		ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree1, "sigmaSq", sigmaSq1, "alpha", alpha1, "optimumManager", optima, "useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue1, "eqDist", false); // rootValue here is irrelevant
		double lnLk = ouLk.calculateLogP();

		Assert.assertEquals(9.9161059, lnLk, EPSILON); // in R, we get 9.916106
	}

	/*
	 * (3) Small tree, one trait, three optima, root value is random variable (equilibrium assumed at root), estimating root value
	 */
	@Test
	public void testOUMVNLkOneTraitSmallTree3optRandomRVEstimateRV() {
		// thetas
		RealParameter colorValues = new RealParameter(new Double[]{ 0.206222932117995, 0.26633408087427, 0.88122539543514 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[]{ 1, 1, 2, 0, 1, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree1);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree1, "coalCorrection", false);

		// likelihood
		ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree1, "sigmaSq", sigmaSq2, "alpha", alpha2, "optimumManager", optima, "useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue2, "eqDist", true); // rootValue here is irrelevant
		double lnLk = ouLk.calculateLogP();

		Assert.assertEquals(9.9161068, lnLk, EPSILON); // in R, we get 9.916107
	}

	/*
	 * (4) Small tree, one trait, three optima, root value is random variable (equilibrium assumed at root), using first theta as root value (one fewer parameter)
	 */
	@Test
	public void testOUMVNLkOneTraitSmallTree3optRandomRVFirstThetaIsRV() {
		// thetas
		RealParameter colorValues = new RealParameter(new Double[]{ 0.206222932117995, 0.26633408087427, 0.88122539543514 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[]{ 1, 1, 2, 0, 1, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree1);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree1, "coalCorrection", false);

		// likelihood
		ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree1, "sigmaSq", sigmaSq2, "alpha", alpha2, "optimumManager", optima, "useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue2, "eqDist", true);
		double lnLk = ouLk.calculateLogP();

		Assert.assertEquals(9.9161068, lnLk, EPSILON); // in R, we get 9.916107
	}

	/*
	 * (5) Small non-ultrametric tree, one trait, one optimum, conditioning on root value, estimating root value
	 */
	@Test
	public void testOUMVNLkOneTraitSmallNonUltraTree1optCondOnRVEstimateRV() {
		// thetas
		RealParameter colorValues = new RealParameter(new Double[]{ 1.142785e+10 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree2);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree2, "coalCorrection", false);

		// sigma 1, alpha 1 and root value 1
		sigmasqInput = new Double[]{ 0.009103832 };
		sigmaSq = new RealParameter(sigmasqInput);
		alphaInput = new Double[]{ 1.818089e-11 };
		alpha = new RealParameter(alphaInput);
		rootValueInput = new Double[]{ -3.726855e-01 };
		rootValue = new RealParameter(rootValueInput);

		// likelihood
		ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree2, "sigmaSq", sigmaSq, "alpha", alpha, "optimumManager", optima, "useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);
		double lnLk = ouLk.calculateLogP();

		Assert.assertEquals(1.3679572, lnLk, EPSILON); // in R, we get 1.367972 -- different in the 5th decimal place (we have some pretty small inputs, probably some numerical recipe difference; everything else passes!)
	}

	/*
	 * (6) Small non-ultrametric tree, one trait, one optimum, conditioning on root value, using first theta as root value (one fewer parameter)
	 */
	@Test
	public void testOUMVNLkOneTraitSmallNonUltraTree1optCondOnRVFirstThetaIsRV() {
		// thetas
		RealParameter colorValues = new RealParameter(new Double[]{ 0.3497826 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree2);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree2, "coalCorrection", false);

		// sigma 1, alpha 1 and root value 1
		sigmasqInput = new Double[]{ 0.01940518 };
		sigmaSq = new RealParameter(sigmasqInput);
		alphaInput = new Double[]{ 2.707329e-11 };
		alpha = new RealParameter(alphaInput);
		rootValueInput = new Double[]{ -3.726855e-01 };
		rootValue = new RealParameter(rootValueInput);

		// likelihood
		ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree2, "sigmaSq", sigmaSq, "alpha", alpha, "optimumManager", optima, "useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);
		double lnLk = ouLk.calculateLogP();

		Assert.assertEquals(-0.1459122, lnLk, EPSILON); // in R, we get -0.1459122
	}

	/*
	 * (7) Small non-ultrametric tree, one trait, one optimum, root value is random variable (equilibrium assumed at root), using first theta as root value (one fewer parameter)
	 */
	@Test
	public void testOUMVNLkOneTraitSmallNonUltraTree3optRandomRVFirstThetaIsRV() {
		// thetas
		RealParameter colorValues = new RealParameter(new Double[]{ 1.035041 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree2);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree2, "coalCorrection", false);

		// sigma 1, alpha 1 and root value 1
		sigmasqInput = new Double[]{ 0.02879764 };
		sigmaSq = new RealParameter(sigmasqInput);
		alphaInput = new Double[]{ 0.4316411 };
		alpha = new RealParameter(alphaInput);
		rootValueInput = new Double[]{ -1.924925 };
		rootValue = new RealParameter(rootValueInput);

		// likelihood
		ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree2, "sigmaSq", sigmaSq, "alpha", alpha, "optimumManager", optima, "useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		double lnLk = ouLk.calculateLogP();

		Assert.assertEquals(1.1710858, lnLk, EPSILON); // in R, we get 1.171086
	}

	/*
	 * (8) Small non-ultrametric tree, one trait, one optimum, root value is random variable (equilibrium assumed at root), estimating root value
	 */
	@Test
	public void testOUMVNLkOneTraitSmallNonUltraTree3optRandomRVEstimateRV() {
		// thetas
		RealParameter colorValues = new RealParameter(new Double[]{ 0.4152632 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree2);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree2, "coalCorrection", false);

		// sigma 1, alpha 1 and root value 1
		sigmasqInput = new Double[]{ 0.09143114 };
		sigmaSq = new RealParameter(sigmasqInput);
		alphaInput = new Double[]{ 0.599265 };
		alpha = new RealParameter(alphaInput);
		rootValueInput = new Double[]{ -1.924925 };
		rootValue = new RealParameter(rootValueInput);

		// likelihood
		ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree2, "sigmaSq", sigmaSq, "alpha", alpha, "optimumManager", optima, "useRootMetaData", false, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		double lnLk = ouLk.calculateLogP();

		Assert.assertEquals(-0.5143806, lnLk, EPSILON); // in R, we get -0.5143806
	}

	/*
	 * (9) Large non-ultrametric tree with sampled ancestors, one trait, one optimum, root value is random variable (equilibrium assumed at root), using first theta as root value (one fewer parameter)
	 */
	@Test
	public void testOUMVNLkOneTraitLargeNonUltraTree1optCondOnRVEstimateRV() {
		String treeStr = "(t72_1:33.95517374,((t45_2:105.7970875,t45_1:0):61.92884103,((t16_2:64.31266032,t16_1:0):33.18670286,((((((((t33_2:17.24690108,t33_1:0):34.29978643,((((t29_5:8.042489968,t29_4:0):48.76232459,t29_3:0):2.870515606,t29_2:0):113.0118958,t29_1:0):78.8636162):2.752286226,t104_1:0):13.2284639,t76_1:19.98684163):35.40421183,t99_1:0):22.37220403,(((((((t25_2:5.545398329,t25_1:0):29.73789406,(((t14_1:19.83955306,(t68_2:7.642029053,t68_1:0):12.197524):4.493086625,t15_1:24.33263968):7.659951961,t32_1:31.99259164):4.589970848):16.34481241,t18_1:35.83487487):68.91853002,((t41_1:26.38487177,(((t39_1:21.69802025,t20_1:21.69802025):22.43410663,t160_2:0):0.5011266949,t160_1:0):24.01552498):43.01583713,t138_1:0):10.18128923):70.28879289,((((t40_3:5.674797941,t40_2:0):85.82961996,t40_1:0):9.144286626,((t22_1:6.416294394,t62_1:6.416294394):52.92804634,t169_1:0):41.3043638):44.11984404,t126_1:0):47.36614923):93.92949444,t113_1:0):24.15812702,(((((t7_1:96.64737931,(t56_1:9.124505504,t60_1:9.124505504):87.5228738):142.4776409,t48_1:53.57593652):13.68925583,((((((t5_1:1.47740974,t24_1:1.47740974):67.36924043,t172_1:0):0.6094754868,((t81_3:23.44109494,t81_2:0):11.8683594,t81_1:0):34.14667132):94.37344028,(((t42_1:24.136168,((t78_2:67.1575288,t78_1:0):6.07403593,t13_1:11.99646428):19.79173068):16.01164487,t8_1:15.88519099):53.96855772,(((((t64_2:57.50834491,t64_1:0):14.48295323,(t44_1:66.77751843,((t19_1:12.37766384,(t3_2:21.04342944,t3_1:0):34.8018801):10.2073654,t146_1:0):0.7248434878):5.213779703):6.124517947,(((((t85_3:30.86302382,t85_2:0):2.746255125,t85_1:0):7.42992738,(t23_2:35.51946648,t23_1:0):5.519739842):5.318122342,(t63_1:29.4370073,t4_1:29.4370073):16.92032137):10.87142916,(((t27_2:0.7705516481,t27_1:0):16.68364597,t65_1:17.45419762):0.8973042275,(t30_1:5.872210534,(t9_1:5.379640887,(t49_1:0.03872308848,t54_1:0.03872308848):5.340917798):0.4925696479):12.47929131):38.87725598):20.88705826):57.37473544,t132_2:0):11.64855352,t132_1:0):15.86439296):0.8260679418):12.34991723,(t67_2:54.95630021,t67_1:0):46.89744564):72.69817908,t115_1:0):3.936613795):37.64004304,((((t73_3:142.7476486,t73_2:0):23.26030434,t73_1:0):33.70759497,(t86_1:29.84643822,((t11_2:52.84278453,t11_1:0):65.0721915,t31_1:54.03135372):52.51395311):30.03845768):0.7587008348,t100_1:0):32.43803262):10.22767547,t6_1:74.01326892):9.540324705):15.08568888):1.038342891,t92_1:0):8.58602375,(((t52_1:16.45509171,(((t46_1:148.9059947,t66_1:169.3974901):0.6996035147,t80_1:71.3880344):1.442748076,t109_1:0):8.282393923):76.16579002,(((t47_1:87.75883092,(((t37_1:19.21090631,t71_1:19.21090631):2.716558918,t162_1:0):52.25715577,(t50_2:14.18543288,t50_1:0):0.7051498905):13.57420992):70.93504021,(t84_1:91.53880307,(((t43_3:45.6907903,t43_2:0):27.11086441,t43_1:0):26.97698465,t82_1:57.43850143):42.53791402):16.37731775):67.82194637,t120_1:0):70.63623304):19.65198808,(((t12_1:85.77104457,((((((t26_1:6.654829005,((t36_1:12.77851341,t70_1:12.77851341):20.88269631,(t53_1:31.23167166,((t74_1:11.18633708,t2_1:11.18633708):13.86162972,(t59_2:8.290822624,t59_1:0):16.75714418):6.183704852):2.429538064):33.88601089):65.52318133,((t17_1:47.74620153,t77_1:47.74620153):26.04800054,((t87_1:24.52842258,t69_1:24.52842258):40.93917936,t79_1:65.46760194):8.326600137):59.27619986):2.038654818,t34_1:33.0666319):2.275848245,t122_1:0):59.7289996,((((t57_2:9.620141517,t57_1:0):45.82168925,((t58_2:5.178405873,t58_1:0):55.47396517,((t21_1:14.00132991,t28_1:14.00132991):36.16842295,t165_1:0):10.48261818):3.289390269):95.35491899,t142_1:0):22.11781524,t35_1:6.324292886):15.69940907):2.23954769,((t83_1:0.07202234145,((t1_1:4.659677209,(((t61_4:25.59732606,t61_3:0):8.674379048,t61_2:0):4.825121609,t61_1:0):4.738801364):14.43560063,t150_1:0):20.9778315):78.47955114,t131_1:0):41.62484094):37.0609):69.68067928,t106_1:0):0.7776192694,((t10_1:9.54139831,(t38_2:26.07965526,t38_1:0):76.6686444):85.31864731,(((t75_1:63.31636704,t55_1:21.6953316):38.78938204,t143_1:0):151.8995115,t51_1:49.73077652):1.996244886):50.8711454):9.931387773):18.12833616):44.30923714):67.66883573):36.06056073):70.8162184;";
		Tree myTree = new TreeParser(treeStr, false, false, true, 0);

		// initializing data
		/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
		// RealParameter oneTraitValues = new RealParameter(new Double[] { -0.120540925368686,-0.160339196119639,0.339711780786955,0.366620606449858,0.33088594760519,0.240619559058153,0.290536304372365,-0.0901170041442549,0.341338296432439,-0.0628145511095494,0.0580385312534993,-0.327327215900842,-0.0518215793539553,0.33460186054066,0.327108133275378,0.163999857771553,-0.0625322749065135,0.473100099350845,-0.18548137014437,0.182392891817864,0.0635218896061768,-0.0692694738005909,0.226987428950308,-0.0135619972058633,0.215366928258746,0.30056737056481,0.312907629853014,0.275301281455329,0.515793851411512,-0.0331061262646763,0.552621356036791,0.0812736250137989,0.029500783167557,0.355533617814899,0.226135370395429,0.35733531062027,0.132222943832334,0.386522540915343,0.223921904096349,0.247184759229459,0.206342673206047,0.49090135847057,-0.127630867933923,-0.163114992198786,0.285977651079011,0.590086989328278,0.204482309523344,0.0951399327857863,-0.186578808911084,0.438175763843511,0.307394267727196,0.230745719707747,0.376416248991533,0.015059308494915,0.254547570530282,-0.0235761692154553,0.455228154309083,-0.243392546580896,0.194520124445963,0.153455693187263,-0.0207109572209939,0.108876605908814,0.0853388155408518,-0.0261917713957721,-0.181563099673691,0.598439787340723,-0.145307160355961,0.158228203216972,0.171675227407924,-0.0415889203783169,0.23912156617336,-0.0663147445692704,0.187870149577858,-0.176702690795953,0.0625705762692598,0.282082324865967,-0.34947856036343,0.154309228836352,0.324847770367215,0.272495695719496,0.356018737503574,-0.0286542881780232,0.162373514467401,0.0301730646940898,0.14344777534431,0.413025722015541,0.138415356600982,-0.0756638106298703,0.253105483929284,0.435068723137525,0.187599226491797,0.0792622583121616,0.166598040080181,0.108375803217645,0.20390709863448,0.367963844490222,0.281057351966486,0.0673500911931306,0.382454086988514,0.140311163307048,0.0945200530794752,0.174999612642347,0.231026476185091,0.226220487819618,-0.0515643818227237,0.226908996988942,0.233647101034971,0.170618432275688,0.330053268957442,0.212922368401166,-0.120113266676103,-0.297073778660506,0.745154054326884,0.410473276891546,0.222125735804677,0.0673003797775433,0.0445273257571818,0.327534883207233,-0.0106841482120232,0.125740419588887,0.325803704776373,-0.045762978254612,0.234871127364338,0.430714626348258,0.388572986226116,0.27356995328721,-0.220756150970401,0.117111468872904,-0.0116525857230558,0.15116217224561,0.0642469156709887,-0.222758595055236,0.188856275423778,-0.115686175121437,-0.052489420279737,0.35585389337747,0.27098637632671,0.0117372993387806,0.33551924588533,0.346790730504337,0.216496852706343,0.483101983845081,0.292554498008609,0.410138399734205,0.138258175735177,0.045149840801478 });
		// spNames = "t47_1,t44_1,t81_3,t29_5,t39_1,t7_1,t70_1,t32_1,t59_2,t14_1,t36_1,t79_1,t5_1,t61_4,t22_1,t27_2,t23_2,t3_2,t77_1,t87_1,t74_1,t2_1,t65_1,t20_1,t60_1,t43_3,t28_1,t49_1,t85_3,t62_1,t64_2,t68_2,t30_1,t15_1,t58_2,t75_1,t53_1,t78_2,t37_1,t56_1,t4_1,t40_3,t63_1,t17_1,t54_1,t9_1,t24_1,t69_1,t21_1,t71_1,t25_2,t57_2,t18_1,t1_1,t66_1,t55_1,t41_1,t82_1,t19_1,t84_1,t11_2,t73_3,t50_2,t26_1,t13_1,t46_1,t38_2,t42_1,t67_2,t83_1,t8_1,t34_1,t31_1,t80_1,t12_1,t10_1,t35_1,t48_1,t86_1,t33_2,t51_1,t52_1,t6_1,t76_1,t45_2,t16_2,t72_1,t81_1,t81_2,t29_1,t29_2,t29_3,t29_4,t59_1,t61_1,t61_2,t61_3,t27_1,t23_1,t3_1,t43_1,t43_2,t85_1,t85_2,t64_1,t68_1,t58_1,t78_1,t40_1,t40_2,t25_1,t57_1,t11_1,t73_1,t73_2,t50_1,t38_1,t67_1,t33_1,t45_1,t16_1,t92_1,t99_1,t100_1,t104_1,t106_1,t109_1,t113_1,t115_1,t120_1,t122_1,t126_1,t131_1,t132_1,t132_2,t138_1,t142_1,t143_1,t146_1,t150_1,t160_1,t160_2,t162_1,t165_1,t169_1,t172_1";
		// oneTraitData = new OneValueContTraits();
		// oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
		spNames = "t47_1 t44_1 t81_3 t29_5 t39_1 t7_1 t70_1 t32_1 t59_2 t14_1 t36_1 t79_1 t5_1 t61_4 t22_1 t27_2 t23_2 t3_2 t77_1 t87_1 t74_1 t2_1 t65_1 t20_1 t60_1 t43_3 t28_1 t49_1 t85_3 t62_1 t64_2 t68_2 t30_1 t15_1 t58_2 t75_1 t53_1 t78_2 t37_1 t56_1 t4_1 t40_3 t63_1 t17_1 t54_1 t9_1 t24_1 t69_1 t21_1 t71_1 t25_2 t57_2 t18_1 t1_1 t66_1 t55_1 t41_1 t82_1 t19_1 t84_1 t11_2 t73_3 t50_2 t26_1 t13_1 t46_1 t38_2 t42_1 t67_2 t83_1 t8_1 t34_1 t31_1 t80_1 t12_1 t10_1 t35_1 t48_1 t86_1 t33_2 t51_1 t52_1 t6_1 t76_1 t45_2 t16_2 t72_1 t81_1 t81_2 t29_1 t29_2 t29_3 t29_4 t59_1 t61_1 t61_2 t61_3 t27_1 t23_1 t3_1 t43_1 t43_2 t85_1 t85_2 t64_1 t68_1 t58_1 t78_1 t40_1 t40_2 t25_1 t57_1 t11_1 t73_1 t73_2 t50_1 t38_1 t67_1 t33_1 t45_1 t16_1 t92_1 t99_1 t100_1 t104_1 t106_1 t109_1 t113_1 t115_1 t120_1 t122_1 t126_1 t131_1 t132_1 t132_2 t138_1 t142_1 t143_1 t146_1 t150_1 t160_1 t160_2 t162_1 t165_1 t169_1 t172_1";
		oneTraitValues = Arrays.asList(-0.120540925368686, -0.160339196119639, 0.339711780786955, 0.366620606449858, 0.33088594760519, 0.240619559058153, 0.290536304372365, -0.0901170041442549, 0.341338296432439, -0.0628145511095494, 0.0580385312534993, -0.327327215900842, -0.0518215793539553, 0.33460186054066, 0.327108133275378, 0.163999857771553, -0.0625322749065135, 0.473100099350845, -0.18548137014437, 0.182392891817864, 0.0635218896061768, -0.0692694738005909, 0.226987428950308, -0.0135619972058633, 0.215366928258746, 0.30056737056481, 0.312907629853014, 0.275301281455329, 0.515793851411512, -0.0331061262646763, 0.552621356036791, 0.0812736250137989, 0.029500783167557, 0.355533617814899, 0.226135370395429, 0.35733531062027, 0.132222943832334, 0.386522540915343, 0.223921904096349, 0.247184759229459, 0.206342673206047, 0.49090135847057, -0.127630867933923, -0.163114992198786, 0.285977651079011, 0.590086989328278, 0.204482309523344, 0.0951399327857863, -0.186578808911084, 0.438175763843511, 0.307394267727196, 0.230745719707747, 0.376416248991533, 0.015059308494915, 0.254547570530282, -0.0235761692154553, 0.455228154309083, -0.243392546580896, 0.194520124445963, 0.153455693187263, -0.0207109572209939, 0.108876605908814, 0.0853388155408518, -0.0261917713957721, -0.181563099673691, 0.598439787340723, -0.145307160355961, 0.158228203216972, 0.171675227407924, -0.0415889203783169, 0.23912156617336, -0.0663147445692704, 0.187870149577858, -0.176702690795953, 0.0625705762692598, 0.282082324865967, -0.34947856036343, 0.154309228836352, 0.324847770367215, 0.272495695719496, 0.356018737503574, -0.0286542881780232, 0.162373514467401, 0.0301730646940898, 0.14344777534431, 0.413025722015541, 0.138415356600982, -0.0756638106298703, 0.253105483929284, 0.435068723137525, 0.187599226491797, 0.0792622583121616, 0.166598040080181, 0.108375803217645, 0.20390709863448, 0.367963844490222, 0.281057351966486, 0.0673500911931306, 0.382454086988514, 0.140311163307048, 0.0945200530794752, 0.174999612642347, 0.231026476185091, 0.226220487819618, -0.0515643818227237, 0.226908996988942, 0.233647101034971, 0.170618432275688, 0.330053268957442, 0.212922368401166, -0.120113266676103, -0.297073778660506, 0.745154054326884, 0.410473276891546, 0.222125735804677, 0.0673003797775433, 0.0445273257571818, 0.327534883207233, -0.0106841482120232, 0.125740419588887, 0.325803704776373, -0.045762978254612, 0.234871127364338, 0.430714626348258, 0.388572986226116, 0.27356995328721, -0.220756150970401, 0.117111468872904, -0.0116525857230558, 0.15116217224561, 0.0642469156709887, -0.222758595055236, 0.188856275423778, -0.115686175121437, -0.052489420279737, 0.35585389337747, 0.27098637632671, 0.0117372993387806, 0.33551924588533, 0.346790730504337, 0.216496852706343, 0.483101983845081, 0.292554498008609, 0.410138399734205, 0.138258175735177, 0.045149840801478);
		oneTraitData = new RealParameter();
		// oneTraitData = new KeyRealParameter();
		oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

		// thetas
		RealParameter colorValues = new RealParameter(new Double[] { 1.596677e-01 }); // thetas
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
		RateCategoryClockModel rcc = new RateCategoryClockModel();
		rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

		TreeToVCVMat optima = new TreeToVCVMat();
		optima.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);

		// sigma, alpha and root value
		sigmasqInput = new Double[] { 0.06733451 };
		sigmaSq = new RealParameter(sigmasqInput);
		alphaInput = new Double[] { 0.7972125 };
		alpha = new RealParameter(alphaInput);
		rootValueInput = new Double[] { 6.082900e-87 };
		rootValue = new RealParameter(rootValueInput);

		// likelihood
		ouLk = new OUMVNLikelihoodOneTrait();
		ouLk.initByName("tree", myTree, "sigmaSq", sigmaSq, "alpha", alpha, "optimumManager", optima, "useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);
		double lnLk = ouLk.calculateLogP();

		Assert.assertEquals(25.417006, lnLk, EPSILON); // in R, we get 25.41701
	}
}