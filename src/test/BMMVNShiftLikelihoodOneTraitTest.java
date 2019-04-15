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

public class BMMVNShiftLikelihoodOneTraitTest {

	double lnLk, lnLk2;
	final static double EPSILON = 1e-6;
	
	/*
	 * Small tree, simple BM. Second likelihood adds root edge. Should match/reduce to BMMVNLikelihoodOneTraitTest
	 */
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.2704762 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0 });
		IntegerParameter rootEdgeColorAssignment = new IntegerParameter(new Integer[] { 0 });
		ColorManager colors = new ColorManager();
		colors.initByName("tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		ColorManager colors2 = new ColorManager();
		colors2.initByName("tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", true, "rootEdgeLength", 1.0, "rootEdgeColorAssignment", rootEdgeColorAssignment);
				
		// initializing data	
		RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9 });
		String spNames = "sp1,sp2,sp3";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		// mean vector
		Double[] meanVectorInput = new Double[] { 4.985714 };
		RealParameter mean = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNShiftLikelihoodOneTrait BMLk = new BMMVNShiftLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "rateManager", colors, "mean", mean, "oneTraitData", oneTraitData);
		lnLk = BMLk.calculateLogP(); // no root edge
		
		BMMVNShiftLikelihoodOneTrait BMLk2 = new BMMVNShiftLikelihoodOneTrait();
		BMLk2.initByName("tree", myTree, "rateManager", colors2, "mean", mean, "oneTraitData", oneTraitData);
		lnLk2 = BMLk2.calculateLogP();
	}

	@Test
	public void testLnLk() {
		Assert.assertEquals(-3.191339, lnLk, EPSILON);
		Assert.assertEquals(-3.577933, lnLk2, EPSILON);
	}
}
