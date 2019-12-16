package test;

import beast.evolution.tree.Node;
import org.apache.commons.math3.linear.*;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import beast.util.TreeParser;
import contraband.PruneUtils;

import java.util.ArrayList;
import java.util.List;

/*
 * In this unit test, we run PruneUtils using a fixed tree with 8 taxa
 * and test the sets methods of vCVMat, aMat, bVec, cMat, dVec, eMat, f.
 * We select two nodes:
 * 1) node number 10, at which there is a 10 * 10 evolutionary rate matrix.
 * 2) node number 3, at which there is a 15 * 15 evolutionary rate matrix.
 */
public class PruneUtilsTest1 {
    final static double EPSILON = 1e-2;
    private RealMatrix vCVMat1; private RealMatrix aMat1; private RealMatrix eMat1; private RealMatrix cMat1; private RealVector dVec1; private double f1;
    private RealMatrix vCVMat2; private RealMatrix aMat2; private RealMatrix eMat2; private RealMatrix cMat2; private RealVector dVec2; private double f2;
    @Before
    public void setUP() throws Exception {
        // tree
        String treeStr = "(((t2:0.3911985469,t3:0.3911985469):0.07316864044,(t5:0.1281555032,t6:0.1281555032):0.3362116841):0.5356328126,((t4:0.1323580755,(t7:0.03532969201,t8:0.03532969201):0.09702838344):0.8127611167,t1:0.9451191922):0.05488080785):0.0;";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);


        // initiate lists
        RealMatrix initRM1 = new Array2DRowRealMatrix(new double[10][10]);
        RealVector iniVec1 = new ArrayRealVector(new double[10]);
        RealMatrix initRM2 = new Array2DRowRealMatrix(new double [15][15]);
        RealVector iniVec2 = new ArrayRealVector(new double[15]);

        List<RealMatrix> vcvMatList = new ArrayList<>(tree.getNodeCount());
        List<RealMatrix> aMatList = new ArrayList<>(tree.getNodeCount());
        List<RealMatrix> eMatList = new ArrayList<>(tree.getNodeCount());
        List<RealMatrix> cMatList = new ArrayList<>(tree.getNodeCount());
        List<RealVector> dVecList = new ArrayList<>(tree.getNodeCount());
        double[] fArrary = new double[tree.getNodeCount()];

        // Node1
        Node node1 = tree.getNode(10);
        for (int i = 0; i < tree.getNodeCount(); i++) {
            vcvMatList.add(i, initRM1);
            aMatList.add(i, initRM1);
            eMatList.add(i, initRM1);
            cMatList.add(i, initRM1);
            dVecList.add(i, iniVec1);
        }

        double[][] evolRateMat2DArray1 = {{9.568333, -0.9333169, 1.452668, -5.078245, 7.790786, 4.1706094, 8.0459809, -3.636380, -5.339318, -5.869372},
                                        {0.000000, 6.7757064, -7.941506, - 9.158809, 3.856068, 0.8813205, 3.8141056, -5.367484, -0.680751, -7.449367},
                                        {0.000000, 0.0000000, 8.998250, -3.441586, 2.810136, 1.8828404, 5.9093484, -7.144000, -4.680547, 5.066157},
                                        {0.000000, 0.0000000, 0.000000, 9.545036, 9.885396, -4.2168053, -9.5077263, -1.709073, 7.156554, 7.900907},
                                        {0.000000, 0.0000000, 0.000000, 0.000000, 6.557058, -7.0577271, -0.4440806, -1.725513, -9.083377, -2.510744},
                                        {0.000000, 0.0000000, 0.000000, 0.000000, 0.000000, 9.6302423, 5.1691908, -2.623091, -1.155999, 3.302304},
                                        {0.000000, 0.0000000, 0.000000, 0.000000, 0.000000, 0.0000000, 2.1640794, -6.951105, 5.978497, -8.103187},
                                        {0.000000, 0.0000000, 0.000000, 0.000000, 0.000000, 0.0000000, 0.0000000,  1.388061, -7.562015, -2.320607},
                                        {0.000000, 0.0000000, 0.000000, 0.000000, 0.000000, 0.0000000, 0.0000000,  0.000000, 5.609480, -4.512327},
                                        {0.000000, 0.0000000, 0.000000, 0.000000, 0.000000, 0.0000000, 0.0000000,  0.000000, 0.000000, 8.146400}};
        RealMatrix evolRateMat1 = new Array2DRowRealMatrix(evolRateMat2DArray1);
        PruneUtils.setVCVMatForBranchBM(node1, evolRateMat1,vcvMatList);
        vCVMat1 = vcvMatList.get(node1.getNr());

        LUDecomposition VMatLUD1 = new LUDecomposition(vCVMat1);
        RealMatrix inverseVCVMat1 = VMatLUD1.getSolver().getInverse();

        PruneUtils.setAMat(node1, inverseVCVMat1, aMatList);
        aMat1 = aMatList.get(node1.getNr());

        RealMatrix phiMat1 = MatrixUtils.createRealIdentityMatrix(10);
        PruneUtils.setEMatOU(node1, inverseVCVMat1, phiMat1.transpose(), eMatList);
        eMat1 = eMatList.get(node1.getNr());

        PruneUtils.setCMatOU(node1, eMat1, phiMat1, cMatList);
        cMat1 = cMatList.get(node1.getNr());

        RealVector omegaVec1 = new ArrayRealVector(new double [10]);
        PruneUtils.setDVecOU(node1, eMat1, omegaVec1, dVecList);
        dVec1 = dVecList.get(node1.getNr());

        double vCVMatDet1 = VMatLUD1.getDeterminant();
        PruneUtils.setF(node1, 10, vCVMatDet1, fArrary);
        f1 = fArrary[node1.getNr()];

        // Node2
        Node node2 = tree.getNode(3);
        for (int i = 0; i < tree.getNodeCount(); i++) {
            vcvMatList.add(i, initRM2);
            aMatList.add(i, initRM2);
            eMatList.add(i, initRM2);
            cMatList.add(i, initRM2);
            dVecList.add(i, iniVec2);
        }
        double[][] evolRateMat2DArray2 = {{8.99825, -5.0782453, -3.441586, 3.856068,  4.1706094, 9.260485, 5.169191, -2.6230910, -1.155999,  3.302304,  5.0895032,  3.361112,  3.1351626,  8.2887637,  2.957870},
                                          {0.00000, 0.4205953, 9.090073, 2.810136,  0.8813205,  8.045981, -5.671841, -6.9511050,  5.978497, -8.103187,  2.5844226, -1.647064, -3.5925352,  2.1746996, -3.603588},
                                          {0.00000, 0.0000000, 8.895393, 9.885396,  1.8828404,  3.814106, -3.636380, -7.2238787, -7.562015, -2.320607,  4.2036480,  5.763917, -6.2461776, -1.7862045, -3.845600},
                                          {0.00000, 0.0000000, 0.000000, 6.557058, -4.2168053,  5.909348, -5.367484, -5.3393180,  1.218960, -4.512327, -9.9875045, -7.942707,  5.6458860, -7.0581062, -5.604647},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  1.4711365, -9.507726, -7.144000, -0.6807510, -5.869372,  6.292801, -0.4936685, -1.302145, -8.1281003,  8.7059961, -2.610223},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  0.0000000,  4.777960, -1.709073, -4.6805472, -7.449367, -1.029673, -5.5976223,  9.699140, -0.6644192, -3.9754220,  9.684384},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  0.0000000,  0.000000,  4.137243,  7.1565543,  5.066157,  6.201287, -2.4036692,  7.861022,  0.2301092, -8.7855886, -6.915954},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  0.0000000,  0.000000,  0.000000,  0.4583117,  7.900907,  6.247790,  2.2554201,  7.729381,  1.9997792,  8.9545388, -8.179120},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  0.0000000,  0.000000,  0.000000,  0.0000000,  3.744628,  5.886846, -2.9640418, -6.498947, -3.3435292,  4.4119255, -7.161862},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  0.0000000,  0.000000,  0.000000,  0.0000000,  0.000000,  4.398317, -7.7772915, -7.386086, -0.2277393, -7.1541141,  3.800142},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  0.0000000,  0.000000,  0.000000,  0.0000000,  0.000000,  0.000000,  2.4361947,  3.062039,  9.0894765,  0.9856931,  2.385130},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  0.0000000,  0.000000,  0.000000,  0.0000000,  0.000000,  0.000000,  0.0000000,  3.435165, -0.3419521,  9.0818248,  7.827882},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  0.0000000,  0.000000,  0.000000,  0.0000000,  0.000000,  0.000000,  0.0000000,  0.000000,  8.9035022,  1.7096671,  3.459982},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  0.0000000,  0.000000,  0.000000,  0.0000000,  0.000000,  0.000000,  0.0000000,  0.000000,  0.0000000,  4.0451028,  4.741555},
                                          {0.00000, 0.0000000, 0.000000, 0.000000,  0.0000000,  0.000000,  0.000000,  0.0000000,  0.000000,  0.000000, 0.0000000,  0.000000,  0.0000000,  0.0000000,  5.211357}};
        RealMatrix evolRateMat2 = new Array2DRowRealMatrix(evolRateMat2DArray2);
        PruneUtils.setVCVMatForBranchBM(node2, evolRateMat2, vcvMatList);
        vCVMat2 = vcvMatList.get(node2.getNr());

        LUDecomposition VMatLUD2 = new LUDecomposition(vCVMat2);
        RealMatrix inverseVCVMat2 = VMatLUD2.getSolver().getInverse();

        PruneUtils.setAMat(node2, inverseVCVMat2, aMatList);
        aMat2 = aMatList.get(node2.getNr());

        RealMatrix phiMat2 = MatrixUtils.createRealIdentityMatrix(15);
        PruneUtils.setEMatOU(node2, inverseVCVMat2, phiMat2.transpose(), eMatList);
        eMat2 = eMatList.get(node2.getNr());

        PruneUtils.setCMatOU(node2, eMat2, phiMat2, cMatList);
        cMat2 = cMatList.get(node2.getNr());

        RealVector omegaVec2 = new ArrayRealVector(new double [15]);
        PruneUtils.setDVecOU(node2, eMat2, omegaVec2, dVecList);
        dVec2 = dVecList.get(node2.getNr());

        double vCVMatDet2 = VMatLUD2.getDeterminant();
        PruneUtils.setF(node2, 15, vCVMatDet2, fArrary);
        f2 = fArrary[node2.getNr()];
    }
    @Test
    public void againstresLikelihood () {
        // test first node
        Assert.assertEquals(-2.276653, vCVMat1.getEntry(5,6), EPSILON);
        Assert.assertEquals(-0.18263775, aMat1.getEntry(4,5), EPSILON);
        Assert.assertEquals(65.1254377, eMat1.getEntry(7,8), EPSILON);
        Assert.assertEquals(0.49655071, cMat1.getEntry(6,4), EPSILON);
        Assert.assertEquals(0.0000, dVec1.getEntry(3), EPSILON);
        Assert.assertEquals(-23.76009, f1, EPSILON);

        // test second node
        Assert.assertEquals(6.679951, vCVMat2.getEntry(5,14), EPSILON);
        Assert.assertEquals(-564.39925, aMat2.getEntry(4,5), EPSILON);
        Assert.assertEquals(-62.343741, eMat2.getEntry(3,2), EPSILON);
        Assert.assertEquals(16675.77, cMat2.getEntry(5,10), EPSILON);
        Assert.assertEquals(0.0000, dVec2.getEntry(3), EPSILON);
        Assert.assertEquals(-16.76498 , f2, EPSILON);
    }
}
