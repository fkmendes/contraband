package test;

import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.Test;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.coalescent.CoalCorrection;

/**
 * @author Fabio K. Mendes
 */

public class CoalCorrectionTest {
	
	final static double EPSILON = 1e-8;
	double[] expected1row, expected2row, expected3row, expected4row, expected5row;
	double[][] correctedPhyloTMat;

	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,(sp4:1.0,sp5:1.0):2.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		String[] spNamesInPhyloTMatOrder = new String[myTree.getLeafNodeCount()];

		// pop sizes
		Double[] popSizesInput = new Double[] { 0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.75, 1.0, 0.5 };
		RealParameter popSizes = new RealParameter(popSizesInput);
		
		CoalCorrection coal = new CoalCorrection();
		coal.initByName("tree", myTree, "popSizes", popSizes);
		correctedPhyloTMat = coal.getCorrectedPhyloTMat(spNamesInPhyloTMatOrder);
		
		//GeneralUtils.display2DArray(correctedPhyloTMat);
		
		expected1row = new double[] { 3.5638913869461915, 2.3059405550002823, 0.8797906714751229, 0.06389138694619145, 0.06389138694619145 };
		expected2row = new double[] { 2.3059405550002823, 3.5638913869461915, 0.8797906714751229, 0.06389138694619145, 0.06389138694619145 };
		expected3row = new double[] { 0.8797906714751229, 0.8797906714751229, 3.5638913869461915, 0.06389138694619145, 0.06389138694619145 };
		expected4row = new double[] { 0.06389138694619145, 0.06389138694619145, 0.06389138694619145, 3.5638913869461915, 1.631559028564498 };
		expected5row = new double[] { 0.06389138694619145, 0.06389138694619145, 0.06389138694619145, 1.631559028564498, 3.5638913869461915 };
	}
	
	@Test
	public void testLnLk() {
		assertArrayEquals(expected1row, correctedPhyloTMat[0], EPSILON);
		assertArrayEquals(expected2row, correctedPhyloTMat[1], EPSILON);
		assertArrayEquals(expected3row, correctedPhyloTMat[2], EPSILON);
		assertArrayEquals(expected4row, correctedPhyloTMat[3], EPSILON);
		assertArrayEquals(expected5row, correctedPhyloTMat[4], EPSILON);
	}
}
