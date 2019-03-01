package test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.evolution.tree.Node;
import beast.util.TreeParser;
import contraband.MVNUtils;

public class MVNGetVcvMatTestDriver {

	public static void main(String[] args) {
		String treeStr = "(((sp39:0.518912972,sp40:0.518912972):19.54206195,((((sp25:3.198513788,(sp32:2.293402763,(sp41:0.1728996412,sp42:0.1728996412):2.120503122):0.9051110254):7.525323533,((sp20:9.577599427,(sp26:2.751623892,(sp30:2.405609293,sp31:2.405609293):0.3460145989):6.825975535):0.05909341512,(sp17:7.221384607,sp18:7.221384607):2.415308235):1.087144479):0.5715875464,sp10:11.29542487):6.453137462,((sp33:2.252609903,sp34:2.252609903):4.989398146,sp16:7.24200805):10.50655428):2.312412597):0.3952852439,(((((sp23:5.605262355,sp24:5.605262355):3.179619681,((sp37:1.329072526,sp38:1.329072526):1.780228265,(sp27:2.543803164,sp28:2.543803164)nd38:0.5654976265):5.675581245):6.165501477,sp6:14.95038351):0.6290423683,((sp11:8.298747349,(sp15:8.099808068,sp12:8.099808068):0.1989392817):3.788483262,(sp21:6.890801228,(sp35:1.989124199,sp36:1.989124199):4.901677029):5.196429383):3.492195269):2.878773551,sp4:18.45819943):1.998060738)0.0;";
		String treeStrNonUltra = "(((sp39:0.8,sp40:0.4):10.0,((((sp25:3.198513788,(sp32:2.293402763,(sp41:0.1728996412,sp42:0.1728996412):2.120503122):0.9051110254):7.525323533,((sp20:9.577599427,(sp26:2.751623892,(sp30:2.405609293,sp31:2.405609293):0.3460145989):6.825975535):0.05909341512,(sp17:7.221384607,sp18:7.221384607):2.415308235):1.087144479):0.5715875464,sp10:11.29542487):6.453137462,((sp33:2.252609903,sp34:2.252609903):4.989398146,sp16:7.24200805):10.50655428):2.312412597):0.3952852439,(((((sp23:5.605262355,sp24:5.605262355):3.179619681,((sp37:1.329072526,sp38:1.329072526):1.780228265,(sp27:2.543803164,sp28:2.543803164)nd38:0.5654976265):5.675581245):6.165501477,sp6:14.95038351):0.6290423683,((sp11:8.298747349,(sp15:8.099808068,sp12:8.099808068):0.1989392817):3.788483262,(sp21:6.890801228,(sp35:1.989124199,sp36:1.989124199):4.901677029):5.196429383):3.492195269):2.878773551,sp4:18.45819943):1.998060738)0.0;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		TreeParser myTreeNonUltra = new TreeParser(treeStrNonUltra, false, false, true, 0);
		int n = myTree.getLeafNodeCount();
		double[][] tMat = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
		
		MVNUtils.populateVcvMatrix(myTree, nodeToRootPaths, tMat, leftLeaves, rightLeaves);
		System.out.println(tMat[7][5]); // cov(sp20, sp17) = 10.8195673283
		
		MVNUtils.populateVcvMatrix(myTreeNonUltra, nodeToRootPaths, tMat, leftLeaves, rightLeaves);
		System.out.println(tMat[24][26]); // cov(sp39, sp40) = 10.3952852439
		System.out.println(tMat[24][24]); // var(sp39) = 11.1952852439
		System.out.println(tMat[26][26]); // var(sp40) = 10.7952852439
		
		for (int i=0; i<tMat.length; i++) {
			System.out.println(Arrays.toString(tMat[i]));
		}
	}

}
