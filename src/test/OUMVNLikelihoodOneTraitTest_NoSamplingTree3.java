package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.GeneralUtils;
import contraband.OUMVNLikelihoodOneTrait_NoTreeSampling;
import contraband.OneValueContTraits;

public class OUMVNLikelihoodOneTraitTest_NoSamplingTree3 {

	double lnLk;
	final static double EPSILON = 1e-5;
	
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((((t35[&Regime=0]:2.336518061,t32[&Regime=0]:2.336518061)[&Regime=0]:28.95257479,t10[&Regime=0]:31.28909285)[&Regime=0]:8.654086516,t18[&Regime=0]:39.94317937)[&Regime=0]:52.28906298,(((t47[&Regime=0]:31.00652286,t9[&Regime=0]:31.00652286)[&Regime=0]:50.20634817,(t43[&Regime=0]:15.06939472,t38[&Regime=0]:15.06939472)[&Regime=0]:66.14347631)[&Regime=0]:10.61662549,(((((t20[&Regime=0]:20.94406932,t14[&Regime=0]:20.94406932)[&Regime=0]:28.09437292,t19[&Regime=0]:49.03844224)[&Regime=0]:31.16991698,(t24[&Regime=0]:54.88723469,(((t50[&Regime=0]:2.534803909,t8[&Regime=0]:2.534803909)[&Regime=0]:35.18774941,t25[&Regime=0]:37.72255332)[&Regime=0]:15.41876911,(t12[&Regime=0]:42.01655137,t5[&Regime=0]:42.01655137)[&Regime=0]:11.12477106)[&Regime=0]:1.745912255)[&Regime=0]:25.32112453)[&Regime=0]:2.610368667,t37[&Regime=0]:82.81872788)[&Regime=0]:6.617999642,(t42[&Regime=0]:81.65977864,(t13[&Regime=0]:5.88018515,t41[&Regime=0]:5.88018515)[&Regime=0]:75.77959349)[&Regime=0]:7.776948892)[&Regime=0]:2.392768999)[&Regime=0]:0.4027458181)[&Regime=0]:7.767757656,((t34[&Regime=0]:80.73867518,((t4[&Regime=0]:14.89974775,t36[&Regime=0]:14.89974775)[&Regime=0]:7.855467399,t7[&Regime=0]:22.75521515)[&Regime=0]:57.98346003)[&Regime=0]:16.48666894,((((((((t29[&Regime=0]:32.9204832,t22[&Regime=0]:32.9204832)[&Regime=0]:13.17504731,t46[&Regime=0]:46.09553051)[&Regime=0]:1.732718052,t40[&Regime=0]:47.82824856)[&Regime=0]:14.51317295,(t28[&Regime=0]:29.85457377,((t33[&Regime=0]:6.373725141,t21[&Regime=0]:6.373725141)[&Regime=0]:1.191235246,t26[&Regime=0]:7.564960387)[&Regime=0]:22.28961339)[&Regime=0]:32.48684774)[&Regime=0]:5.177695495,t48[&Regime=0]:67.51911701)[&Regime=0]:2.445324178,(t39[&Regime=0]:56.9237382,((t2[&Regime=0]:5.876590264,t44[&Regime=0]:5.876590264)[&Regime=0]:19.06403767,t23[&Regime=0]:24.94062793)[&Regime=0]:31.98311027)[&Regime=0]:13.04070299)[&Regime=0]:0.3095854321,(((t11[&Regime=0]:13.30542076,t49[&Regime=0]:13.30542076)[&Regime=0]:14.69428372,t45[&Regime=0]:27.99970449)[&Regime=0]:1.437902517,t31[&Regime=0]:29.43760701)[&Regime=0]:40.83641961)[&Regime=0]:11.48412211,((((t16[&Regime=0]:30.59346099,(t30[&Regime=0]:0.03406798076,t1[&Regime=0]:0.03406798076)[&Regime=0]:30.55939301)[&Regime=0]:21.47527084,(t17[&Regime=0]:50.41024027,t15[&Regime=0]:50.41024027)[&Regime=0]:1.658491566)[&Regime=0]:14.63237622,(t3[&Regime=0]:10.35007739,t27[&Regime=0]:10.35007739)[&Regime=0]:56.35103066)[&Regime=0]:2.944577857,t6[&Regime=0]:69.64568591)[&Regime=0]:12.11246283)[&Regime=0]:15.46719539)[&Regime=0]:2.774655878)[&Regime=0]:0;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data		
		RealParameter oneTraitValues = new RealParameter(new Double[] { 1.46036403223613,1.46684738870839,1.51293402720201,1.38645983622075,1.52360274878184,1.38604516251543,1.34674375692537,1.57061239621065,1.48962576974379,1.48653875766911,1.47882999523348,1.42987407092074,1.500830609038,1.52573554902885,1.35101467964535,1.3754167761214,1.39992186810313,1.47268752575019,1.30540867869118,1.43286309430355,1.35747497708841,1.42191922014755,1.34837424935406,1.5102597527159,1.4259500689268,1.43069041000211,1.52900698556248,1.39163255825086,1.35189645355176,1.40593882022866,1.40399457029073,1.46245982107446,1.42363339595256,1.35978645751607,1.43131557265059,1.4768092041642,1.36373666377089,1.50038165340804,1.48327796797891,1.52949573694115,1.42292947537042,1.4666710783077,1.39114425864009,1.42708261181083,1.38061320641582,1.50472419662598,1.35762598538961,1.47786657493344,1.40219215464719,1.37165305130316 });
		String spNames = "t37,t42,t24,t19,t12,t5,t18,t25,t10,t47,t9,t20,t14,t43,t38,t13,t41,t50,t8,t35,t32,t34,t6,t48,t39,t17,t15,t40,t46,t29,t22,t16,t28,t31,t45,t23,t7,t4,t36,t11,t49,t3,t27,t26,t33,t21,t2,t44,t30,t1";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues, "spNames", spNames);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.01507664 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
		
		// alpha
		Double[] alphaInput = new Double[] { 1.988807 };
		RealParameter alpha = new RealParameter(alphaInput);
		
		// root value
		Double[] rootValueInput = new Double[] { 6.082900e-87 };
		RealParameter rootValue = new RealParameter(rootValueInput);
		
		// theta
		Double[] thetaInput = { 1.435158 };
		RealParameter theta = new RealParameter(thetaInput);
		
		// likelihood 1 (condition on rootValue, theta_0 as parameter)
		OUMVNLikelihoodOneTrait_NoTreeSampling OULk = new OUMVNLikelihoodOneTrait_NoTreeSampling();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "theta", theta, "nOptima", 1,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", false);
		lnLk = OULk.calculateLogP();
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(69.16449, lnLk, EPSILON); 
	}
}
