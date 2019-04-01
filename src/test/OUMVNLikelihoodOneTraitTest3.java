package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.GeneralUtils;
import contraband.OUMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

public class OUMVNLikelihoodOneTraitTest3 {

	double lnLk;
	final static double EPSILON = 1e-5;
	
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((((t35[&Regime=0]:2.336518061,t32[&Regime=0]:2.336518061)[&Regime=0]:28.95257479,t10[&Regime=0]:31.28909285)[&Regime=0]:8.654086516,t18[&Regime=0]:39.94317937)[&Regime=0]:52.28906298,(((t47[&Regime=0]:31.00652286,t9[&Regime=0]:31.00652286)[&Regime=0]:50.20634817,(t43[&Regime=0]:15.06939472,t38[&Regime=0]:15.06939472)[&Regime=0]:66.14347631)[&Regime=0]:10.61662549,(((((t20[&Regime=0]:20.94406932,t14[&Regime=0]:20.94406932)[&Regime=0]:28.09437292,t19[&Regime=0]:49.03844224)[&Regime=0]:31.16991698,(t24[&Regime=0]:54.88723469,(((t50[&Regime=0]:2.534803909,t8[&Regime=0]:2.534803909)[&Regime=0]:35.18774941,t25[&Regime=0]:37.72255332)[&Regime=0]:15.41876911,(t12[&Regime=0]:42.01655137,t5[&Regime=0]:42.01655137)[&Regime=0]:11.12477106)[&Regime=0]:1.745912255)[&Regime=0]:25.32112453)[&Regime=0]:2.610368667,t37[&Regime=0]:82.81872788)[&Regime=0]:6.617999642,(t42[&Regime=0]:81.65977864,(t13[&Regime=0]:5.88018515,t41[&Regime=0]:5.88018515)[&Regime=0]:75.77959349)[&Regime=0]:7.776948892)[&Regime=0]:2.392768999)[&Regime=0]:0.4027458181)[&Regime=0]:7.767757656,((t34[&Regime=0]:80.73867518,((t4[&Regime=0]:14.89974775,t36[&Regime=0]:14.89974775)[&Regime=0]:7.855467399,t7[&Regime=0]:22.75521515)[&Regime=0]:57.98346003)[&Regime=0]:16.48666894,((((((((t29[&Regime=0]:32.9204832,t22[&Regime=0]:32.9204832)[&Regime=0]:13.17504731,t46[&Regime=0]:46.09553051)[&Regime=0]:1.732718052,t40[&Regime=0]:47.82824856)[&Regime=0]:14.51317295,(t28[&Regime=0]:29.85457377,((t33[&Regime=0]:6.373725141,t21[&Regime=0]:6.373725141)[&Regime=0]:1.191235246,t26[&Regime=0]:7.564960387)[&Regime=0]:22.28961339)[&Regime=0]:32.48684774)[&Regime=0]:5.177695495,t48[&Regime=0]:67.51911701)[&Regime=0]:2.445324178,(t39[&Regime=0]:56.9237382,((t2[&Regime=0]:5.876590264,t44[&Regime=0]:5.876590264)[&Regime=0]:19.06403767,t23[&Regime=0]:24.94062793)[&Regime=0]:31.98311027)[&Regime=0]:13.04070299)[&Regime=0]:0.3095854321,(((t11[&Regime=0]:13.30542076,t49[&Regime=0]:13.30542076)[&Regime=0]:14.69428372,t45[&Regime=0]:27.99970449)[&Regime=0]:1.437902517,t31[&Regime=0]:29.43760701)[&Regime=0]:40.83641961)[&Regime=0]:11.48412211,((((t16[&Regime=0]:30.59346099,(t30[&Regime=0]:0.03406798076,t1[&Regime=0]:0.03406798076)[&Regime=0]:30.55939301)[&Regime=0]:21.47527084,(t17[&Regime=0]:50.41024027,t15[&Regime=0]:50.41024027)[&Regime=0]:1.658491566)[&Regime=0]:14.63237622,(t3[&Regime=0]:10.35007739,t27[&Regime=0]:10.35007739)[&Regime=0]:56.35103066)[&Regime=0]:2.944577857,t6[&Regime=0]:69.64568591)[&Regime=0]:12.11246283)[&Regime=0]:15.46719539)[&Regime=0]:2.774655878)[&Regime=0]:0;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data		
		String oneTraitValues = "t37=0.499504982506037,t42=0.493600773470607,t24=0.38929832284563,t19=0.325244394827358,t12=0.373350590070522,t5=0.41319496893665,t18=0.441849920184112,t25=0.529666176159169,t10=0.479036987313547,t47=0.481990640613883,t9=0.48927832812467,t20=0.434436964548657,t14=0.379435120689744,t43=0.389724893953509,t38=0.362931638929068,t13=0.3666929875179,t41=0.353607267558753,t50=0.492116715241339,t8=0.371243960115714,t35=0.478121086631717,t32=0.578741325960089,t34=0.430496718172667,t6=0.431313899417636,t48=0.4909254403674,t39=0.441193352930789,t17=0.397518292080972,t15=0.359299983058545,t40=0.447477695803669,t46=0.417142120463927,t29=0.464908396116619,t22=0.458476344489581,t16=0.517695584505652,t28=0.477023183396581,t31=0.508275599364981,t45=0.446099967281363,t23=0.489595426373694,t7=0.427623992504684,t4=0.433084084276035,t36=0.446391776928927,t11=0.377326313531013,t49=0.371770458754527,t3=0.425097509765054,t27=0.360342879490491,t26=0.447615875101796,t33=0.48856947342418,t21=0.390604842461911,t2=0.445069161958563,t44=0.352221522955928,t30=0.496528933914092,t1=0.468867719555597";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.01225534 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
		
		// alpha
		Double[] alphaInput = new Double[] { 1.996261 };
		RealParameter alpha = new RealParameter(alphaInput);
		
		// root value
		Double[] rootValueInput = new Double[] { 8.764411e-88 };
		RealParameter rootValue = new RealParameter(rootValueInput);
		
		// theta
		Double[] thetaInput = { 4.357571e-01 };
		RealParameter theta = new RealParameter(thetaInput);
		
		// likelihood 1 (condition on rootValue, theta_0 as parameter)
		OUMVNLikelihoodOneTrait OULk = new OUMVNLikelihoodOneTrait();
		OULk.initByName("tree", myTree, "sigmasq", sigmasq, "alpha", alpha, "theta", theta, "nOptima", 1,
				"useRootMetaData", true, "oneTraitData", oneTraitData, "rootValue", rootValue, "eqDist", true);
		lnLk = OULk.calculateLogP();
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(74.42928, lnLk, EPSILON); 
	}
}
