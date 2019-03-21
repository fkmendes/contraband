package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.BMMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

/*
 * Same as Test2, but tree was made non-ultrametric, negative mean
 */
public class BMMVNLikelihoodOneTraitTest3 {

	double lnLk;
	final static double EPSILON = 1e-6;
	
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((((t36:0.1,t15:0.1):0.1,((t32:0.1,t38:0.1):0.1,t30:0.1):0.1):0.1,t48:0.1):0.1,(((((t29:0.1,((t13:0.1,t10:0.1):0.1,((t49:0.1,t44:0.1):0.1,(t12:0.1,(t42:0.1,t28:0.1):0.1):0.1):0.1):0.1):0.1,((t14:0.1,t35:1):0.1,t45:0.1):0.1):0.1,(t8:0.1,(t6:0.1,t47:0.1):0.1):0.1):0.1,((t19:0.1,((t50:0.1,t1:0.1):0.1,t41:0.1):0.1):0.1,t25:0.1):0.1):0.1,t7:0.1):0.1):0.1,((((t40:5.547816064,t39:5.547816064):69.34457458,(t5:63.2303403,(t3:8.905705974,(t2:8.507831451,t24:8.507831451):0.3978745225):54.32463433):11.66205034):2.166829643,t43:77.05922028):20.83631922,(((t21:35.66443979,(t31:19.79691334,t22:19.79691334):15.86752645):43.35619033,(t16:0.798839978,t27:0.798839978):78.22179014):7.146714057,(((t17:11.8481771,(t4:8.654399431,(t33:4.946918061,t20:4.946918061):3.70748137):3.193777673):60.81563855,(t37:26.39721604,(t26:19.7151034,t9:19.7151034):6.682112638):46.26659961):2.831035132,((t23:4.660327926,t18:4.660327926):67.53247256,((t46:6.150543246,t11:6.150543246):29.58028484,t34:35.73082809):36.4619724):3.302050294):10.6724934):11.72819533):2.104460492):0;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data
		String oneTraitValues = "t36=-0.612880234099481,t15=-0.774278307760373,t32=-0.365516580183378,t38=-0.667458066688047,t30=-0.888717015322935,t48=-0.803604790253531,t29=-0.0563311378876053,t13=0.0649577490538714,t10=-0.213010012835546,t49=0.152854042801119,t44=-0.132222808057845,t12=0.0235699353780309,t42=-0.439443395335116,t28=-0.523115202769663,t14=-0.113413125024881,t35=0.358035099753072,t45=-0.00522561294425657,t8=-0.0492414003489426,t6=-0.290791396373298,t47=-0.0296273685981681,t19=0.182016005682877,t50=0.01766515271966,t1=-0.103410675035095,t41=0.201895614166552,t25=-0.0537343372936915,t7=-0.117563112293728,t40=1.56849130745049,t39=2.73215430891677,t5=-0.538035355088687,t3=-6.31428328740763,t2=-5.69844534887218,t24=-5.9785671742528,t43=-4.64885312184872,t21=4.15367600184535,t31=1.86064361298819,t22=-1.39508204240147,t16=-1.47535379622077,t27=-1.22789051901583,t17=9.08846944938777,t4=7.11666675658611,t33=8.76387924245766,t20=8.89863489298473,t37=1.80636056927312,t26=5.61238617067368,t9=7.03481064594507,t23=3.0247015107113,t18=1.14134052131152,t46=4.21879153861552,t11=4.70952287305256,t34=3.88910154150986";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "traitValues", oneTraitValues);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 0.188522 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
		
		// mean vector
//		Double[] meanVectorInput = new Double[] { 0.059292 };
		Double[] meanVectorInput = new Double[] { -0.415788 };
		RealParameter mean = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNLikelihoodOneTrait BMLk = new BMMVNLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "sigmasq", sigmasq, "mean", mean, "oneTraitData", oneTraitData);
		lnLk = BMLk.calculateLogP();
		System.out.println(lnLk);
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-54.634967, lnLk, EPSILON); 
	}
}
