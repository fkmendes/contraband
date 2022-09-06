package testdrivers;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import java.util.Arrays;
import java.util.List;

public class BMPruneMissingDataTestDriver {
    public static void main(String[] args) {
        // tree
        String treeStr = "((((t5:0.1347986037,t6:0.1347986037):1.147066172,t1:1.281864775):0.1391117368,(t2:0.880414795,(t7:0.1265605926,t8:0.1265605926):0.7538542024):0.5405617172):0.4254821739,((t9:0.06373797515,t10:0.06373797515):0.4034950816,(t3:0.3802203837,t4:0.3802203837):0.08701267311):1.379225629);";
        //String treeStr = "((((t5:0.1347986037,t6:0.1347986037):1.147066172,t1:1.281864775):0.1391117368,(t2:0.880414795,(t7:0.1265605926,t8:0.1265605926):0.7538542024):0.5405617172):0.4254821739,(t3:0.3802203837,t4:0.3802203837):1.379225629);";
        Tree tree = new TreeParser(treeStr, false, false, true, 0);

        int nSpecies = tree.getLeafNodeCount();
        boolean[] isSpeciesToIgnore = new boolean[2 * nSpecies - 1];
        //isSpeciesToIgnore[2] = true;
        //isSpeciesToIgnore[3] = true;
        isSpeciesToIgnore[3] = true;
        isSpeciesToIgnore[4] = true;
        isSpeciesToIgnore[1] = true;
        isSpeciesToIgnore[9] = true;
        //isSpeciesToIgnore[8] = true;

        boolean[] nodeHasMissingData = new boolean[2 * nSpecies - 1];
        int[] speciesToIgnoreIndex = new int[2 * nSpecies - 1];

        prune(tree.getRoot(), isSpeciesToIgnore, nodeHasMissingData, speciesToIgnoreIndex);

        System.out.println( Arrays.toString(isSpeciesToIgnore));

        System.out.println( Arrays.toString(nodeHasMissingData));

        System.out.println( Arrays.toString(speciesToIgnoreIndex));

    }

    public static void prune(Node node, boolean[] isSpeciesToIgnore, boolean[] nodeHasMissingData, int[] speciesToIgnoreIndex){
        List<Node> children = node.getChildren();
        for (Node child : children) {
            if(!child.isLeaf()){
                prune(child, isSpeciesToIgnore, nodeHasMissingData, speciesToIgnoreIndex);
                Node c1 = child.getChild(0);
                Node c2 = child.getChild(1);
                if(isSpeciesToIgnore[c1.getNr()] && isSpeciesToIgnore[c2.getNr()]){
                   isSpeciesToIgnore[child.getNr()] = true;
                   nodeHasMissingData[node.getNr()] = true;

                    Node sib = node.getChild(0);
                    int sibNr = sib.getNr();
                    if(sibNr  == child.getNr()){
                        speciesToIgnoreIndex[node.getNr()] = node.getChild(1).getNr();
                    } else {
                        speciesToIgnoreIndex[node.getNr()] = sibNr;
                    }
                }

                if(isSpeciesToIgnore[c1.getNr()] && !isSpeciesToIgnore[c2.getNr()]){
                    nodeHasMissingData[child.getNr()] = true;
                    speciesToIgnoreIndex[child.getNr()] = c2.getNr();
                }

                if(!isSpeciesToIgnore[c1.getNr()] && isSpeciesToIgnore[c2.getNr()]){
                    nodeHasMissingData[child.getNr()] = true;
                    speciesToIgnoreIndex[child.getNr()] = c1.getNr();
                }


            }
        }

    }

}
