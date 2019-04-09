package test;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.ColorManager;

public class ColorManagerTest {

	double[][] colorValuesMat;
	
	double[] expected1row, expected2row;
	
	@Before
	public void setUp() throws Exception {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):2.0,(sp4:2.5,sp5:2.5):1.5);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
						
		// initializing data
		RealParameter colorValues = new RealParameter(new Double[] { 0.2, 0.4, 0.6, 0.8, 1.0 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 3, 3, 0, 0, 0, 2, 1, 4, 0 });
				
		ColorManager colors = new ColorManager();
		colors.initByName("tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments);
		
		colorValuesMat = colors.getSpColorValuesMatrix();
		
		expected1row = new double[] {
				((0.8*1.0) + (0.6*1.0) + (0.4*2.0)),
				(0.6*1.0) + (0.4*2.0),
				(0.4*2.0),
				0.0, 0.0 };
		expected2row = new double[] { 0.0, 0.0, 0.0, 0.0 };
	}

	@Test
	public void test() {
		assertArrayEquals(expected1row, colorValuesMat[0], 1E-6);
	}

}
