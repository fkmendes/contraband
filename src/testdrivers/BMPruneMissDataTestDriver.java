package testdrivers;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

import java.util.Arrays;
import java.util.List;

public class BMPruneMissDataTestDriver {

    public static void main(String[] args) {
        // tree
        String treeStr = "(((t4:4.157956267,(t5:3.645264322,t6:3.645264322):0.5126919446):3.853008658,(t2:5.632320957,t3:5.632320957):2.378643968):4.138629419,(((t8:2.935100956,(t9:1.420145405,t10:1.420145405):1.514955551):0.5094521248,t7:3.444553081):7.252722915,t1:10.697276):1.452318347);";
        Tree tree = new TreeParser(treeStr, false, false, true, 0);
        boolean[] hasMissingDataSpecies = new boolean[tree.getNodeCount()];
        int[] missingDataSpeciesIndex = new int[tree.getNodeCount()];
        boolean[] speciesToIgnore = new boolean[]{true, true, false, true, false, true, false, false, false, false};

        prune(tree.getRoot(), speciesToIgnore, missingDataSpeciesIndex, hasMissingDataSpecies);
        System.out.println(Arrays.toString(hasMissingDataSpecies));

        System.out.println(Arrays.toString(missingDataSpeciesIndex));
    }

    private static void prune(Node node, boolean[] speciesToIgnore, int[] missingDataSpeciesIndex, boolean[] hasMissingDataSpecies){
        List<Node> children = node.getChildren();
        for (Node child : children) {
            boolean leftChild =false ;
            boolean rightChild = false;
            if(child.isLeaf()){
                if(speciesToIgnore[child.getNr()]){
                    leftChild = true;
                } else {
                    missingDataSpeciesIndex[node.getNr()] = child.getNr();
                }
            }
            else{
                prune(child, speciesToIgnore, missingDataSpeciesIndex, hasMissingDataSpecies);
                if(hasMissingDataSpecies[child.getNr()]){
                    rightChild = true;
                }
            }

            if( leftChild|| rightChild){
                hasMissingDataSpecies[node.getNr()] = true;
            }

        }
    }
}
