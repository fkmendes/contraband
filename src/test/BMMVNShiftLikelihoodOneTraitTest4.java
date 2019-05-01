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

public class BMMVNShiftLikelihoodOneTraitTest4 {

	double lnLk;
	final static double EPSILON = 1e-4;
	
	/*
	 * Large tree, with fossils, one-rate BM (on shift-BM class)
	 */
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.1160941, 0.01914707 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 1, 1, 0, 0, 0 });
		ColorManager colors = new ColorManager();
		colors.initByName("tree", myTree, "nTraits", 1, "nColors", 2, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
				
		// initializing data	
		RealParameter oneTraitValues = new RealParameter(new Double[] { -0.9291812, -0.7312343, -1.6712572 });
		String spNames = "sp1,sp2,sp3";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// mean vector
		Double[] meanVectorInput = new Double[] { -1.125558 };
		RealParameter mean = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "mean", mean, "oneTraitData", oneTraitData);
		lnLk = BMLk.calculateLogP();
		System.out.println(lnLk); // -0.8583677116744495
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-0.8583676, lnLk, EPSILON); 
	}
}
