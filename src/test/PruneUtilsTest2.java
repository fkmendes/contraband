package test;

import beast.evolution.tree.Node;
import org.apache.commons.math3.linear.*;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import beast.util.TreeParser;
import contraband.prunelikelihood.PruneUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Rong Zhang
 */

/*
 * In this unit test, we run PruneUtils using a fixed tree with 8 taxa
 * and test the sets methods of vCVMat, aMat, bVec, cMat, dVec, eMat, f.
 * We select two nodes:
 * 1) node number 10, at which there is a 10 * 10 evolutionary rate matrix.
 * 2) node number 3, at which there is a 15 * 15 evolutionary rate matrix.
 */
public class PruneUtilsTest2 {
    final static double EPSILON = 1e-2;
    RealMatrix lMat; RealVector mVec; double r;
    @Before
    public void setUP() throws Exception {
        // tree
        String treeStr = "((t3:0.2403694672,t4:0.2403694672):1.759630533,(t1:1.819708683,t2:1.819708683):0.1802913169):0.0;";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);


        // initiate lists of aMat, cMat, eMat, aPlusLMat, invAlusLMat,
        // lMat, mVec and traitValues
        int nTraits = 4;
        int nodeCount = tree.getNodeCount();
        List<RealMatrix> aMatList = new ArrayList<>(nodeCount);
        List<RealMatrix> aPlusLListList = new ArrayList<>(nodeCount);
        List<RealMatrix> invAPlusLList = new ArrayList<>(nodeCount);
        List<RealMatrix> lMatList = new ArrayList<>(nodeCount);
        List<RealVector> mVecList = new ArrayList<>(nodeCount);
        List<RealVector> traitValuesVec = new ArrayList<>(nodeCount) ;

        RealMatrix iniRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        RealVector iniVec = new ArrayRealVector(new double [nTraits]);
        for (int i = 0; i < nodeCount; i ++) {
            aMatList.add(i, iniRM);
            mVecList.add(i, iniVec);
            aPlusLListList.add(i,iniRM);
            invAPlusLList.add(i,iniRM);
            lMatList.add(i,iniRM);
            traitValuesVec.add(i, iniVec);
        }
        // initialize arrays of f and r
        double [] fArr = new double[nodeCount];
        double [] rArr = new double[nodeCount];

        // set trait values
        traitValuesVec.set(0, new ArrayRealVector(new double [] {10.577067, 3.268922, -4.028851, 8.857763}));
        traitValuesVec.set(1, new ArrayRealVector(new double [] {-2.627771, 5.437578, -2.347264, 8.229367}));
        traitValuesVec.set(2, new ArrayRealVector(new double [] {24.2162869, 0.1494845, 9.4296224, 4.6112014}));
        traitValuesVec.set(3, new ArrayRealVector(new double [] {18.113242, 4.154542, 6.458512, 2.495724}));

        // set aMat list
        aMatList.set(0, new Array2DRowRealMatrix(new double [][] {
                {-3.106567e-03, -0.005346505,  0.006544133,  2.667865e-05},
                {-5.346505e-03, -0.019053565,  0.013482222, -2.936384e-03},
                {6.544133e-03,  0.013482222, -0.027464139,  8.788890e-03},
                {2.667865e-05, -0.002936384,  0.008788890, -1.435137e-02}
        }));
        aMatList.set(1, new Array2DRowRealMatrix(new double [][] {
                {-3.106567e-03, -0.005346505,  0.006544133,  2.667865e-05},
                {-5.346505e-03, -0.019053565,  0.013482222, -2.936384e-03},
                {6.544133e-03,  0.013482222, -0.027464139,  8.788890e-03},
                {2.667865e-05, -0.002936384,  0.008788890, -1.435137e-02}
        }));
        aMatList.set(2, new Array2DRowRealMatrix(new double [][] {
                {-0.0235181588, -0.04047553,  0.04954213,  0.0002019698},
                {-0.0404755287, -0.14424435,  0.10206670, -0.0222297942},
                {0.0495421279,  0.10206670, -0.20791631,  0.0665359828},
                {0.0002019698, -0.02222979,  0.06653598, -0.1086465677}
        }));
        aMatList.set(3, new Array2DRowRealMatrix(new double [][] {
                {-0.0235181588, -0.04047553,  0.04954213,  0.0002019698},
                {-0.0404755287, -0.14424435,  0.10206670, -0.0222297942},
                {0.0495421279,  0.10206670, -0.20791631,  0.0665359828},
                {0.0002019698, -0.02222979,  0.06653598, -0.1086465677}
        }));

        aMatList.set(4, new Array2DRowRealMatrix(new double [][] {
                {-3.212633e-03, -0.005529048,  0.006767565,  2.758953e-05},
                {-5.529048e-03, -0.019704101,  0.013942539, -3.036640e-03},
                {6.767565e-03,  0.013942539, -0.028401833,  9.088964e-03},
                {2.758953e-05, -0.003036640,  0.009088964, -1.484136e-02}
        }));
        aMatList.set(5, new Array2DRowRealMatrix(new double [][] {
                {-0.0313550724, -0.05396312, 0.06605096, 0.0002692719},
                {-0.0539631161, -0.19231064, 0.13607820, -0.0296373884},
                {0.0660509618, 0.13607820, -0.27719988, 0.0887076483},
                {0.0002692719, -0.02963739, 0.08870765, -0.1448506675}
        }));

        // set cMat list
        List<RealMatrix> cMatList = aMatList;

        // set eMat list
        List<RealMatrix> eMatList = new ArrayList<>(nodeCount);
        eMatList.add(0, new Array2DRowRealMatrix(new double [][] {
                {6.213134e-03, 0.010693010, -0.01308827, -5.335731e-05},
                {1.069301e-02,  0.038107131, -0.02696444,  5.872768e-03},
                {-1.308827e-02, -0.026964445,  0.05492828, -1.757778e-02},
                {-5.335731e-05,  0.005872768, -0.01757778,  2.870275e-02}
        }));
        eMatList.add(1, new Array2DRowRealMatrix(new double [][] {
                {6.213134e-03,  0.010693010, -0.01308827, -5.335731e-05},
                {1.069301e-02,  0.038107131, -0.02696444,  5.872768e-03},
                {-1.308827e-02, -0.026964445,  0.05492828, -1.757778e-02},
                {-5.335731e-05,  0.005872768, -0.01757778,  2.870275e-02}
        }));
        eMatList.add(2, new Array2DRowRealMatrix(new double [][] {
                {0.0470363176,  0.08095106, -0.09908426, -0.0004039396},
                {0.0809510574,  0.28848871, -0.20413339,  0.0444595884},
                {-0.0990842558, -0.20413339,  0.41583261, -0.1330719655},
                {-0.0004039396,  0.04445959, -0.13307197,  0.2172931354}
        }));
        eMatList.add(3, new Array2DRowRealMatrix(new double [][] {
                {0.0470363176,  0.08095106, -0.09908426, -0.0004039396},
                {0.0809510574,  0.28848871, -0.20413339,  0.0444595884},
                {-0.0990842558, -0.20413339,  0.41583261, -0.1330719655},
                {-0.0004039396,  0.04445959, -0.13307197,  0.2172931354}
        }));

        eMatList.add(4, new Array2DRowRealMatrix(new double [][] {
                {6.425266e-03,  0.011058096, -0.01353513, -5.517906e-05},
                {1.105810e-02,  0.039408203, -0.02788508,  6.073279e-03},
                {-1.353513e-02, -0.027885078,  0.05680367, -1.817793e-02},
                {-5.517906e-05,  0.006073279, -0.01817793,  2.968273e-02}
        }));
        eMatList.add(5, new Array2DRowRealMatrix(new double [][] {
                {0.0627101449, 0.10792623, -0.1321019, -0.0005385437},
                {0.1079262322, 0.38462128, -0.2721564, 0.0592747769},
                {-0.1321019236, -0.27215639, 0.5543998, -0.1774152965},
                {-0.0005385437, 0.05927478, -0.1774153, 0.2897013349}
        }));

        // set f array
        fArr[0] = -7.993676; fArr[1] = -7.993676;
        fArr[2] = -12.042185; fArr[3] = -12.042185;
        fArr[4] = -11.975040; fArr[5] = -7.418469;

        List<Node> internalNodes = tree.getInternalNodes();
        for (Node aNode : internalNodes) {
            PruneUtils.setLMatForIntNode(aNode, aMatList, cMatList, eMatList, lMatList, aPlusLListList, invAPlusLList);
            PruneUtils.setMVecForIntNodeBM(aNode, invAPlusLList, eMatList, mVecList, traitValuesVec);
            PruneUtils.setRForIntNode(aNode, aMatList, aPlusLListList, invAPlusLList, mVecList, traitValuesVec, rArr, fArr);
        }

        lMat = lMatList.get(6);
//        System.out.println("printing lMat at root:");
//        GeneralUtils.displayRealMatrix(lMat);
//        System.out.println("printing mVec at root:");
        mVec = mVecList.get(6);
//        GeneralUtils.displayRealVector(mVec);
//        System.out.println("printing r at root:");
        r = rArr[6];
//        System.out.println(r);
    }
    @Test
    public void againstresLikelihood () {
        Assert.assertEquals(0.03555612, lMat.getEntry(1,2), EPSILON);
        Assert.assertEquals(0.5198989, mVec.getEntry(3), EPSILON);
        Assert.assertEquals(-53.72213, r, EPSILON);
    }
}
