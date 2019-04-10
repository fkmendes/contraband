package test;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import contraband.ColorManager;

public class ColorManagerTest {

	public static double[][] colorValuesMat;
	public static double[] expected1row, expected2row, expected3row, expected4row, expected5row;
	public static double sp1, sp3, sp12, expectedSp1, expectedSp3, expectedSp12;
	
	@BeforeClass
	public static void setUp() throws Exception {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):2.0,(sp4:2.5,sp5:2.5):1.5);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
								
		// initializing data
		RealParameter colorValues = new RealParameter(new Double[] { 0.2, 0.4, 0.6, 0.8, 1.0 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 3, 4, 0, 0, 0, 2, 1, 4, 0 });
						
		ColorManager colors = new ColorManager();
		colors.initByName("tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments);
				
		colorValuesMat = colors.getSpColorValuesMat();
		
		expected1row = new double[] {
				((0.8*1.0) + (0.6*1.0) + (0.4*2.0)),
				((0.6*1.0) + (0.4*2.0)),
				0.4*2.0,
				0.0, 0.0 };
		expected2row = new double[] {
				((0.6*1.0) + (0.4*2.0)),
				((1.0*1.0) + (0.6*1.0) + (0.4*2.0)),
				(0.4*2.0),
				0.0, 0.0 };
		expected3row = new double[] {
				0.4*2.0,
				0.4*2.0,
				((0.2*2.0) + (0.4*2.0)),
				0.0, 0.0 };
		expected4row = new double[] {
				0.0, 0.0, 0.0,
				((0.2*2.5) + (1.5*1.0)),
				1.5*1.0 };
		expected5row = new double[] {
				0.0, 0.0, 0.0,
				1.5*1.0,
				((0.2*2.5) + (1.5*1.0))};
		
		Node sp1Node = myTree.getNode(0);
		Node sp3Node = myTree.getNode(2);
		Node sp12Node = myTree.getNode(5);
		
		sp1 = colors.getNodeColorValue(sp1Node);
		sp3 = colors.getNodeColorValue(sp3Node);
		sp12 = colors.getNodeColorValue(sp12Node);
	}
	
	@Test
	public void test1() {
		assertArrayEquals(expected1row, colorValuesMat[0], 1E-10);
		assertArrayEquals(expected2row, colorValuesMat[1], 1E-10);
		assertArrayEquals(expected3row, colorValuesMat[2], 1E-10);
		assertArrayEquals(expected4row, colorValuesMat[3], 1E-10);
		assertArrayEquals(expected5row, colorValuesMat[4], 1E-10);
	}

	@Test
	public void test2() {
		assertEquals(0.8, sp1, 0.0);
		assertEquals(0.2, sp3, 0.0);
		assertEquals(0.6, sp12, 0.0);
	}
}
