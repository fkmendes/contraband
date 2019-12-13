package testdrivers;

import beast.util.TreeParser;
import contraband.GeneralUtils;
import contraband.PruneUtils;
import org.apache.commons.math3.linear.*;
import beast.evolution.tree.Node;
import java.util.ArrayList;
import java.util.List;

public class PruneUtilsTestDriver2 {
    public static void main(String[] args) {
        // initialize bifurcating tree
        // set branch lengths to 1.0 for simplicity
        // tree
        String treeStr =  "((t3:1.0,t4:1.0):1.0,(t1:1.0,t2:1.0):1.0):0.0;";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        // initiate matrix lists
        int nTraits = 4;
        int nodeCount = tree.getNodeCount();
        List<RealMatrix> aMatList = new ArrayList<>(nodeCount);
        RealMatrix iniRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        List<RealVector> mVecList = new ArrayList<>(nodeCount);
        RealVector iniVec = new ArrayRealVector(new double [nTraits]);
        for (int i = 0; i < nodeCount; i ++) {
            aMatList.add(i, iniRM);
            mVecList.add(i, iniVec);
        }
        List<RealMatrix> cMatList = aMatList;
        List<RealMatrix> eMatList = aMatList;
        List<RealMatrix> lMatList = aMatList;
        List<RealMatrix> aPlusLListList = aMatList;
        List<RealMatrix> invAPlusLList = aMatList;
        double [] fArr = new double[nodeCount];
        double [] rArr = fArr;

        // the node to calculate
        Node aIntNode = tree.getNode(4);


        System.out.println("Printing lMat at tip:");
        RealMatrix cMat = new Array2DRowRealMatrix(new double [][]
                {{1.0, 2.0},{3.0, 4.0}}
        );
        RealMatrix lMat = PruneUtils.getLMatForLeaf(cMat);
        GeneralUtils.displayRealMatrix(lMat);

        System.out.println("Printing lMat at this internal node:");
        PruneUtils.setLMatForIntNode(aIntNode, aMatList, cMatList, eMatList, lMatList, aPlusLListList, invAPlusLList);
        GeneralUtils.displayRealMatrix(lMatList.get(aIntNode.getNr()));


        System.out.println("Printing mVec at tip:");
        RealMatrix eMat = new Array2DRowRealMatrix(new double [][]
                {{1.0, 2.0},{3.0, 4.0}}
        );
        RealVector traitValues = new ArrayRealVector(new double []
                {0.8, 0.5, 0.4,0.2});
        RealVector mVec = PruneUtils.getMVecForLeafBM(eMat, traitValues);
        GeneralUtils.displayRealVector(mVec);

        System.out.println("Printing mVec at this internal node:");
        PruneUtils.setMVecForIntNodeBM(aIntNode, invAPlusLList, eMatList, mVecList);
        GeneralUtils.displayRealVector(mVecList.get(aIntNode.getNr()));

        RealMatrix aMat = new Array2DRowRealMatrix(new double [][]
                {{1.0, 2.0},{3.0, 4.0}}
        );
        double f = 1.2;
        double r = PruneUtils.getRForLeafBM(aMat,traitValues, traitValues, f);
        System.out.println("Printing r at leaf:" + r);


        PruneUtils.setRForIntNode(aIntNode, aPlusLListList, invAPlusLList, mVecList, rArr, fArr);
        System.out.println("Printing r at this internal node:" + rArr[aIntNode.getNr()]);



    }
}
