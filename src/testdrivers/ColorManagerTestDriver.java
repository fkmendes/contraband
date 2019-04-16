package testdrivers;

import java.util.Arrays;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import contraband.ColorManager;

public class ColorManagerTestDriver {

	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):2.0,(sp4:2.5,sp5:2.5):1.5);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
				
		// initializing data
		RealParameter colorValues = new RealParameter(new Double[] { 0.2, 0.4, 0.6, 0.8, 1.0 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 3, 3, 0, 0, 0, 2, 1, 4, 0 });
		
		ColorManager colors = new ColorManager();
		colors.initByName("nTraits", 1, "maxNColors", 5, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		double[][] colorValuesMat = colors.getSpColorValuesMatOneTrait();
		for (int i=0; i<colorValuesMat.length; ++i) {
			System.out.println(Arrays.toString(colorValuesMat[i]));
		}
		
		Node sp1Node = myTree.getNode(0);
		System.out.println(colors.getNodeColorValue(sp1Node, 0));
		
		Node sp3Node = myTree.getNode(2);
		System.out.println(colors.getNodeColorValue(sp3Node, 0));
	}

}
