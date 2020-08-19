package testdrivers;

import beast.util.TreeParser;
import contraband.utils.GeneralUtils;
import contraband.prunelikelihood.PruneUtils;
import org.apache.commons.math3.linear.*;

import java.util.ArrayList;
import java.util.List;

/*
 * This test driver runs PruneUtils class using a tree with four taxa,
 * and returns results of sets methods for lMat, mVec and r
 * at both a tip and an internal node.
 */

public class PruneUtilsTestDriver2 {
    public static void main(String[] args) {
        // initialize bifurcating tree
        String treeStr =  "(t1:0.6740016011,(t2:0.6049435608,t3:0.6049435608):0.0690580403):0.0;";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        // initiate matrix lists
        int nTraits = 4;
        int nodeCount = tree.getNodeCount();
        List<RealMatrix> aMatList = new ArrayList<>(nodeCount);
        List<RealMatrix> aPlusLListList = new ArrayList<>(nodeCount);
        List<RealMatrix> invAPlusLList = new ArrayList<>(nodeCount);
        List<RealMatrix> lMatList = new ArrayList<>(nodeCount);
        RealMatrix iniRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        List<RealVector> mVecList = new ArrayList<>(nodeCount);
        RealVector iniVec = new ArrayRealVector(new double [nTraits]);
        List<RealVector> traitValuesVec = new ArrayList<>(nodeCount) ;
        for (int i = 0; i < nodeCount; i ++) {
            aMatList.add(i, iniRM);
            mVecList.add(i, iniVec);
            aPlusLListList.add(i,iniRM);
            invAPlusLList.add(i,iniRM);
            lMatList.add(i,iniRM);
            traitValuesVec.add(i, iniVec);
        }

        double [] fArr = new double[nodeCount];
        double [] rArr = new double[nodeCount];

        // set trait values
        traitValuesVec.set(0, new ArrayRealVector(new double [] {15.508678, 2.783347, 5.825099, 6.283068}));
        traitValuesVec.set(1, new ArrayRealVector(new double [] {-14.932802, 4.835027, -7.601838, 3.817582}));
        traitValuesVec.set(2, new ArrayRealVector(new double [] {-24.6147856, 11.1887244, -12.3152627, 0.4615493}));

        // set aMat list
        aMatList.set(0, new Array2DRowRealMatrix(new double [][] {
                {-8.387291e-03, -0.014434804,  0.01766823,  7.202858e-05},
                {-1.443480e-02, -0.051441923,  0.03640009, -7.927821e-03},
                {1.766823e-02,  0.036400088, -0.07414928,  2.372875e-02},
                {7.202858e-05, -0.007927821,  0.02372875, -3.874667e-02}
        }));
        aMatList.set(1, new Array2DRowRealMatrix(new double [][] {
                {-9.344752e-03, -0.01608263,  0.01968517,  8.025109e-05},
                {-1.608263e-02, -0.05731434,  0.04055538, -8.832830e-03},
                {1.968517e-02,  0.04055538, -0.08261388,  2.643754e-02},
                {8.025109e-05, -0.00883283,  0.02643754, -4.316984e-02}
        }));
        aMatList.set(2, new Array2DRowRealMatrix(new double [][] {
                {-9.344752e-03, -0.01608263,  0.01968517,  8.025109e-05},
                {-1.608263e-02, -0.05731434,  0.04055538, -8.832830e-03},
                {1.968517e-02,  0.04055538, -0.08261388,  2.643754e-02},
                {8.025109e-05, -0.00883283,  0.02643754, -4.316984e-02}
        }));
        aMatList.set(3, new Array2DRowRealMatrix(new double [][] {
                {-0.0818593646, -0.14088267,  0.1724407,  0.0007029938},
                {-0.1408826723, -0.50206954,  0.3552623, -0.0773749699},
                {0.1724406721,  0.35526229, -0.7236917,  0.2315909727},
                { 0.0007029938, -0.07737497,  0.2315910, -0.3781647651}
        }));

        // set cMat list
        List<RealMatrix> cMatList = aMatList;

        // set eMat list
        List<RealMatrix> eMatList = new ArrayList<>(nodeCount);
        eMatList.add(0, new Array2DRowRealMatrix(new double [][] {
                {0.0167745812,  0.02886961, -0.03533646, -0.0001440572},
                {0.0288696088,  0.10288385, -0.07280018,  0.0158556412},
                {-0.0353364587, -0.07280018,  0.14829855, -0.0474575096},
                {-0.0001440572,  0.01585564, -0.04745751,  0.07749333992}
        }));
        eMatList.add(1, new Array2DRowRealMatrix(new double [][] {
                {0.0186895032,  0.03216525, -0.03937033, -0.0001605022},
                {0.0321652527,  0.11462867, -0.08111076,  0.0176656605},
                {-0.0393703336, -0.08111076,  0.16522775, -0.0528750771},
                {-0.0001605022,  0.01766566, -0.05287508,  0.0863396828}
        }));
        eMatList.add(2, new Array2DRowRealMatrix(new double [][] {
                {0.0186895032,  0.03216525, -0.03937033, -0.0001605022},
                {0.0321652527,  0.11462867, -0.08111076,  0.0176656605},
                {-0.0393703336, -0.08111076,  0.16522775, -0.0528750771},
                {-0.0001605022,  0.01766566, -0.05287508,  0.0863396828}
        }));
        eMatList.add(3, new Array2DRowRealMatrix(new double [][] {
                {0.163718729,  0.2817653, -0.3448813, -0.001405988},
                {0.281765345,  1.0041391, -0.7105246,  0.154749940},
                {-0.344881344, -0.7105246,  1.4473834, -0.463181945},
                {-0.001405988,  0.1547499, -0.4631819,  0.756329530}
        }));

        // set f array
        fArr[0] = -10.055786; fArr[1] = -9.839592; fArr[2] = -9.839592; fArr[3] = -5.499216;


        // block for matrix L
        System.out.println("Printing lMat at tip 0:");
        GeneralUtils.displayRealMatrix(PruneUtils.getLMatForLeaf(cMatList.get(0)));
        System.out.println("Printing lMat at tip 1:");
        GeneralUtils.displayRealMatrix(PruneUtils.getLMatForLeaf(cMatList.get(1)));
        System.out.println("Printing lMat at tip 2:");
        GeneralUtils.displayRealMatrix(PruneUtils.getLMatForLeaf(cMatList.get(2)));

        System.out.println("Printing lMat at this internal node 3:");
        PruneUtils.setLMatForIntNode(tree.getNode(3), aMatList, cMatList, eMatList, lMatList, aPlusLListList, invAPlusLList);
        GeneralUtils.displayRealMatrix(lMatList.get(3));

        System.out.println("Printing lMat at this internal node 4:");
        PruneUtils.setLMatForIntNode(tree.getNode(4), aMatList, cMatList, eMatList, lMatList, aPlusLListList, invAPlusLList);
        GeneralUtils.displayRealMatrix(lMatList.get(4));

        // block for vector m
        System.out.println("Printing mVec at tip 0:");
        GeneralUtils.displayRealVector(PruneUtils.getMVecForLeafBM(eMatList.get(0), traitValuesVec.get(0)));
        System.out.println("Printing mVec at tip 1:");
        GeneralUtils.displayRealVector(PruneUtils.getMVecForLeafBM(eMatList.get(1), traitValuesVec.get(1)));
        System.out.println("Printing mVec at tip 2:");
        GeneralUtils.displayRealVector(PruneUtils.getMVecForLeafBM(eMatList.get(2), traitValuesVec.get(2)));

        System.out.println("Printing mVec at this internal node 3:");
        PruneUtils.setMVecForIntNodeBM(tree.getNode(3), invAPlusLList, eMatList, mVecList, traitValuesVec);
        GeneralUtils.displayRealVector(mVecList.get(3));
        System.out.println("Printing mVec at this internal node 4:");
        PruneUtils.setMVecForIntNodeBM(tree.getNode(4), invAPlusLList, eMatList, mVecList, traitValuesVec);
        GeneralUtils.displayRealVector(mVecList.get(4));

        // block for double r
        System.out.println("Printing r at tip 0:" + PruneUtils.getRForLeafBM(aMatList.get(0),traitValuesVec.get(0), fArr[0]));
        System.out.println("Printing r at tip 1:" + PruneUtils.getRForLeafBM(aMatList.get(1),traitValuesVec.get(1), fArr[1]));
        System.out.println("Printing r at tip 2:" + PruneUtils.getRForLeafBM(aMatList.get(2),traitValuesVec.get(2), fArr[2]));

        PruneUtils.setRForIntNode(tree.getNode(3), aMatList, aPlusLListList, invAPlusLList, mVecList, traitValuesVec, rArr, fArr);
        System.out.println("Printing r at this internal node 3:" + rArr[3]);
        PruneUtils.setRForIntNode(tree.getNode(4), aMatList, aPlusLListList, invAPlusLList, mVecList, traitValuesVec, rArr, fArr);
        System.out.println("Printing r at this internal node 4:" + rArr[4]);

    }
}
