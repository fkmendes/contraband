package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.evolution.tree.Node;
import beast.util.TreeParser;

import java.util.ArrayList;
import java.util.List;
import contraband.MVNUtils;

public class MVNPhyloTMatTest {
	
	final static double EPSILON = 1e-4;
	private double[][] phyloTMat, phyloTMat2;
	
	@Before
	public void setUp() throws Exception {
		
		// Test 1: Ultrametric tree
		String treeStr = "(((sp39:0.518912972,sp40:0.518912972):19.54206195,((((sp25:3.198513788,(sp32:2.293402763,(sp41:0.1728996412,sp42:0.1728996412):2.120503122):0.9051110254):7.525323533,((sp20:9.577599427,(sp26:2.751623892,(sp30:2.405609293,sp31:2.405609293):0.3460145989):6.825975535):0.05909341512,(sp17:7.221384607,sp18:7.221384607):2.415308235):1.087144479):0.5715875464,sp10:11.29542487):6.453137462,((sp33:2.252609903,sp34:2.252609903):4.989398146,sp16:7.24200805):10.50655428):2.312412597):0.3952852439,(((((sp23:5.605262355,sp24:5.605262355):3.179619681,((sp37:1.329072526,sp38:1.329072526):1.780228265,(sp27:2.543803164,sp28:2.543803164)nd38:0.5654976265):5.675581245):6.165501477,sp6:14.95038351):0.6290423683,((sp11:8.298747349,(sp15:8.099808068,sp12:8.099808068):0.1989392817):3.788483262,(sp21:6.890801228,(sp35:1.989124199,sp36:1.989124199):4.901677029):5.196429383):3.492195269):2.878773551,sp4:18.45819943):1.998060738)0.0;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		int n = myTree.getLeafNodeCount();
		phyloTMat = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		String[] spNamesInPhyloTMatOrder = new String[myTree.getLeafNodeCount()];
		
		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, phyloTMat, leftLeaves, rightLeaves, spNamesInPhyloTMatOrder);
		
		//System.out.println("cov(sp20, sp17) = " + tMat[7][5]); 		// cov(sp20, sp17) = 10.8195673283
		//System.out.println("cov(sp39, sp40) = " + tMat[24][26]); 	// cov(sp24, sp26) = 19.9373471939
		//System.out.println("cov(sp39, sp39) = " + tMat[24][24]); 	// cov(sp24, sp24) = 20.456260165899998
		//System.out.println("cov(sp40, sp40) = " + tMat[26][26]); 	// cov(sp26, sp26) = 20.456260165899998
		
		// Test 2: non-Ultrametric tree
		String treeStrNonUltra = "(((sp39:0.8,sp40:0.4):10.0,((((sp25:3.198513788,(sp32:2.293402763,(sp41:0.1728996412,sp42:0.1728996412):2.120503122):0.9051110254):7.525323533,((sp20:9.577599427,(sp26:2.751623892,(sp30:2.405609293,sp31:2.405609293):0.3460145989):6.825975535):0.05909341512,(sp17:7.221384607,sp18:7.221384607):2.415308235):1.087144479):0.5715875464,sp10:11.29542487):6.453137462,((sp33:2.252609903,sp34:2.252609903):4.989398146,sp16:7.24200805):10.50655428):2.312412597):0.3952852439,(((((sp23:5.605262355,sp24:5.605262355):3.179619681,((sp37:1.329072526,sp38:1.329072526):1.780228265,(sp27:2.543803164,sp28:2.543803164)nd38:0.5654976265):5.675581245):6.165501477,sp6:14.95038351):0.6290423683,((sp11:8.298747349,(sp15:8.099808068,sp12:8.099808068):0.1989392817):3.788483262,(sp21:6.890801228,(sp35:1.989124199,sp36:1.989124199):4.901677029):5.196429383):3.492195269):2.878773551,sp4:18.45819943):1.998060738)0.0;";
		TreeParser myTreeNonUltra = new TreeParser(treeStrNonUltra, false, false, true, 0);
		int n2 = myTreeNonUltra.getLeafNodeCount();
		phyloTMat2 = new double[n2][n2];
		double[] nodeToRootPaths2 = new double[myTreeNonUltra.getNodeCount()];
		List<Node> leftLeaves2 = new ArrayList<>();
		List<Node> rightLeaves2 = new ArrayList<>();
		String[] spOrderInTMat2 = new String[myTree.getLeafNodeCount()];
		
		MVNUtils.populateTMatrix(myTreeNonUltra, nodeToRootPaths2, phyloTMat2, leftLeaves2, rightLeaves2, spOrderInTMat2);
		
		// System.out.println(spOrderInTMat); // the order they appear is the order in the newick string
		// System.out.println("cov(sp20, sp17) = " + tMat2[7][5]); 	// cov(sp20, sp17) = 10.8195673283
		// System.out.println("cov(sp39, sp40) = " + tMat2[24][26]); 	// cov(sp24, sp26) = 10.3952852439
		// System.out.println("cov(sp39, sp39) = " + tMat2[24][24]); 	// cov(sp24, sp24) = 11.1952852439
		// System.out.println("cov(sp40, sp40) = " + tMat2[26][26]); 	// cov(sp26, sp26) = 10.7952852439
		
	}
	
	@Test // Test1
	public void againstUltPhyloTMat() {
		Assert.assertEquals(10.81957, phyloTMat[7][5], EPSILON); 
		Assert.assertEquals(19.93735, phyloTMat[24][26], EPSILON);
		Assert.assertEquals(20.45626, phyloTMat[24][24], EPSILON);
		Assert.assertEquals(20.45626, phyloTMat[26][26], EPSILON);
	}
	
	@Test // Test2
	public void againstNonUltPhyloTMat() {
		Assert.assertEquals(10.81957, phyloTMat2[7][5], EPSILON); 
		Assert.assertEquals(10.39529, phyloTMat2[24][26], EPSILON);
		Assert.assertEquals(11.19529, phyloTMat2[24][24], EPSILON);
		Assert.assertEquals(10.79529, phyloTMat2[26][26], EPSILON);
	}
	
}
