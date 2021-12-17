package testdrivers;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

import java.util.Arrays;
import java.util.List;

public class BMPruneMissDataTestDriver {

    public static void main(String[] args) {
        // define a complete tree
        String treeStr = "(((t4:4.157956267,(t5:3.645264322,t6:3.645264322):0.5126919446):3.853008658,(t2:5.632320957,t3:5.632320957):2.378643968):4.138629419,(((t8:2.935100956,(t9:1.420145405,t10:1.420145405):1.514955551):0.5094521248,t7:3.444553081):7.252722915,t1:10.697276):1.452318347);";
        Tree tree = new TreeParser(treeStr, false, false, true, 0);

        // some species has missing data denoted by true.
        boolean[] speciesToIgnore = new boolean[]{false, false, true, false, true, true, false, false, true, false};

        // arrays to populate:
        // (1) if an internal node has missing data descendants
        // (2) if yes, the node index of the species at the tip will be recorded, so that the internal node takes the trait values at that tip node
        boolean[] hasMissingDataSpecies = new boolean[tree.getNodeCount()];
        int[] missingDataSpeciesIndex = new int[tree.getNodeCount()];

        // calculate the branch length in time as a scalar of the variance under BM model
        double[] nodeVariance = new double[tree.getNodeCount()];

        // using a recursive algorithm to populate the arrays and calculate the variance
        prune(tree.getRoot(), speciesToIgnore, missingDataSpeciesIndex, hasMissingDataSpecies, nodeVariance);

        // print out the results and compare to the expected
        // the tips are all FALSE
        // the TRUE nodes include 10 11 12 15
        System.out.println("hasMissingDataSpecies " + "\n" +
                "t1: " + hasMissingDataSpecies[0] + "\n" +
                "t10: " + hasMissingDataSpecies[1] + "\n" +
                "t2: " + hasMissingDataSpecies[2] + "\n" +
                "t3: " + hasMissingDataSpecies[3] + "\n" +
                "t4: " + hasMissingDataSpecies[4] + "\n" +
                "t5: " + hasMissingDataSpecies[5] + "\n" +
                "t6: " + hasMissingDataSpecies[6] + "\n" +
                "t7: " + hasMissingDataSpecies[7] + "\n" +
                "t8: " + hasMissingDataSpecies[8] + "\n" +
                "t9: " + hasMissingDataSpecies[9] + "\n" +
                "internal node 10: " + hasMissingDataSpecies[10] + "\n" +
                "internal node 11: " + hasMissingDataSpecies[11] + "\n" +
                "internal node 12: " + hasMissingDataSpecies[12] + "\n" +
                "internal node 13: " + hasMissingDataSpecies[13] + "\n" +
                "internal node 14: " + hasMissingDataSpecies[14] + "\n" +
                "internal node 15: " + hasMissingDataSpecies[15] + "\n" +
                "internal node 16: " + hasMissingDataSpecies[16] + "\n" +
                "internal node 17: " + hasMissingDataSpecies[17] + "\n"
        );

        System.out.println("missingDataSpeciesIndex " + "\n" +
                "internal node 10: " + missingDataSpeciesIndex[10] + "\n" +
                "internal node 11: " + missingDataSpeciesIndex[11] + "\n" +
                "internal node 12: " + missingDataSpeciesIndex[12] + "\n" +
                "internal node 13: " + missingDataSpeciesIndex[13] + "\n" +
                "internal node 14: " + missingDataSpeciesIndex[14] + "\n" +
                "internal node 15: " + missingDataSpeciesIndex[15] + "\n" +
                "internal node 16: " + missingDataSpeciesIndex[16] + "\n" +
                "internal node 17: " + missingDataSpeciesIndex[17] + "\n");

        System.out.println("nodeVariance in branch length (time) " + "\n" +
                "t1: " + nodeVariance[0] + "\t" + "Expect= 10.697276" + "\n" +
                "t10: " + nodeVariance[1] + "\t" + "Expect= 1.4201454049999995" + "\n" +
                "t2: " + nodeVariance[2] + "\t" + "Expect= 0.0" + "\n" +
                "t3: " + nodeVariance[3] + "\t" + "Expect= 5.632320957" + "\n" +
                "t4: " + nodeVariance[4] + "\t" + "Expect= 0.0" + "\n" +
                "t5: " + nodeVariance[5] + "\t" + "Expect= 0.0" + "\n" +
                "t6: " + nodeVariance[6] + "\t" + "Expect= 3.645264322" + "\n" +
                "t7: " + nodeVariance[7] + "\t" + "Expect= 3.4445530810000005" + "\n" +
                "t8: " + nodeVariance[8] + "\t" + "Expect= 0.0" + "\n" +
                "t9: " + nodeVariance[9] + "\t" + "Expect= 1.4201454049999995" + "\n" +
                "internal node 10: " + nodeVariance[10] + "\t" + "Expect= 4.157956267" + "\n" +
                "internal node 11: " + nodeVariance[11] + "\t" + "Expect= 8.010964925" + "\n" +
                "internal node 12: " + nodeVariance[12] + "\t" + "Expect= 8.010964925" + "\n" +
                "internal node 13: " + nodeVariance[13] + "\t" + "Expect= 8.1441118815" + "\n" +
                "internal node 14: " + nodeVariance[14] + "\t" + "Expect= 2.2250282535" + "\n" +
                "internal node 15: " + nodeVariance[15] + "\t" + "Expect= 2.7344803785000007" + "\n" +
                "internal node 16: " + nodeVariance[16] + "\t" + "Expect= 8.777081518078093" + "\n" +
                "internal node 17: " + nodeVariance[17] + "\t" + "Expect= 6.2735743697581965" + "\n" +
                "root: " + nodeVariance[18] + "\t" + "Expect= 3.5437511035976987"
        );

        // backup method for calculating node variance
        //populateNodeVariance(tree.getRoot(), speciesToIgnore, missingDataSpeciesIndex, hasMissingDataSpecies, nodeVariance);
        //System.out.println(Arrays.toString(nodeVariance));


        /*
        t1 = 10.697276 <- parent height <- tip node
        t10 = 1.4201454049999995 <- parent height <- tip node
        t2 = 0.0 <- missing data species
        t3 = 5.632320957 <- parent height <- tip node itself, its sibling is a missing data species
        t4 = 0.0 <- missing data species
        t5 = 0.0 <- missing data species
        t6 = 3.645264322 <- tip node itself, its sibling is a missing data species
        t7 = 3.4445530810000005 <- parent height <- tip node
        t8 = 0.0 <- missing data species
        t9 = 1.4201454049999995 <- parent height <- tip node
        node 10 = (4.157956267 - 3.645264322) + t6 <-  branch length itself + its child  <- an internal node with missing data species
        node 11 = 8.010964925 = (8.010964925 - 4.157956267) + node 10  <-  branch length itself + its child  <- an internal node with missing data species
        node 12 = 8.010964925 = (8.010964925 - 5.632320957) + t3  <-  branch length itself + its child  <- an internal node with missing data species
        node 13 = 4.138629419 + ((node 11 * node 12)/ (node 11 + node 12)) <- a regular internal node
        node 14 = (2.935100956 - 1.420145405) + (1.420145405 * 1.420145405 / (1.420145405 + 1.420145405)) <- a regular internal node
        node 15 = (3.444553081 - 2.935100956) + node 14 <- branch length itself + its child  <- an internal node with missing data species
        node 16 = (10.697276 - 3.444553081) + (node 15 * node 7)/(node 15 + node 7) <- a regular internal node
        node 17 =  1.452318347 + ((node 16 + t1)/(node 16 * t1)) <- a regular internal node
        node 18 = 3.5437511035976987 = (node 13 * node 17)/(node 13 + node 17) <- root calculation
         */

    }

    private static void prune(Node node, boolean[] speciesToIgnore, int[] missingDataSpeciesIndex, boolean[] hasMissingDataSpecies, double[] nodeVariance){
        // initialization
        boolean nodeToIgnore = false;

        // use the node number as index
        int nodeIndex = node.getNr();

        // variance for this node
        nodeVariance[nodeIndex] = node.getLength();

        // get the child nodes of this node
        List<Node> children = node.getChildren();

        // iterate child nodes
        for (Node child : children) {

            int childIndex = child.getNr();

            // CASE1: child is a tip
            if (child.isLeaf()) {
                // if this node has a direct descant that is a species with missing data
                // the flag of this node will be set to TRUE
                nodeToIgnore = speciesToIgnore[childIndex] | nodeToIgnore;

                if (!speciesToIgnore[childIndex]) {
                    // if this child node has trait values
                    // the branch length will be calculated and assigned to this child node
                    nodeVariance[childIndex] = child.getLength();
                }
            }
            // CASE2: child is an internal node
            else {
                // first calculate the length of branch above this child
                nodeVariance[childIndex] = child.getLength();
                // second prune the subtree that is rooted at the child until tips are reached
                prune(child, speciesToIgnore, missingDataSpeciesIndex, hasMissingDataSpecies, nodeVariance);
            }
            // if this node has missing data species
            // we will record the node index of the descant that has trait values
            if(nodeToIgnore){
                missingDataSpeciesIndex[nodeIndex] = childIndex;
            }
        }

        // set the flag
        hasMissingDataSpecies[nodeIndex]  = nodeToIgnore;

        if (nodeToIgnore) {
            // combine the variance
            double compoundVariance = nodeVariance[nodeIndex] + nodeVariance[missingDataSpeciesIndex[nodeIndex]];
            nodeVariance[nodeIndex] = compoundVariance;
        } else {
            // maximum likelihood estimate
            double vc1 = nodeVariance[node.getChild(0).getNr()];
            double vc2 = nodeVariance[node.getChild(1).getNr()];
            nodeVariance[nodeIndex] = nodeVariance[nodeIndex] + ((vc1 * vc2) / (vc1 + vc2));
        }
    }


    private static void populateNodeVariance(Node node, boolean[] speciesToIgnore, int[] missingDataSpeciesIndex, boolean[] hasMissingDataSpecies, double[] nodeVariance){
        List<Node> children = node.getChildren();
        int nodeIndex = node.getNr();

        nodeVariance[nodeIndex] = node.getLength();

        for (Node child : children) {
            int childIndex = child.getNr();

            if(child.isLeaf()) {
                if (!speciesToIgnore[childIndex]) {
                    nodeVariance[childIndex] = child.getLength();
                }
            } else {
                nodeVariance[childIndex] = child.getLength();

                populateNodeVariance(child, speciesToIgnore, missingDataSpeciesIndex, hasMissingDataSpecies, nodeVariance);

            }
        }
        if (hasMissingDataSpecies[nodeIndex]) {
            double compoundVariance = nodeVariance[nodeIndex] + nodeVariance[missingDataSpeciesIndex[nodeIndex]];
            nodeVariance[nodeIndex] = compoundVariance;

        } else {
            double vc1 = nodeVariance[node.getChild(0).getNr()];
            double vc2 = nodeVariance[node.getChild(1).getNr()];

            nodeVariance[nodeIndex] = nodeVariance[nodeIndex] + ((vc1 * vc2) / (vc1 + vc2));
        }
    }
}
