package test;

import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.prunelikelihood.OUPruneLikelihood;
import contraband.prunelikelihood.OUPruneUtils;
import org.apache.commons.math3.linear.*;
import org.junit.Assert;
import org.junit.Test;

/*
 * This class tests methods in OUPruneUtils.class
 */
public class OUPruneUtilsTest {
    final static double EPSILON = 1e-7;

    String treeStr = "((A:3.0058179,B:3.0058179):4.350951,C:7.3567689);";
    Tree tree = new TreeParser(treeStr, false, false, true, 0);

    int nTraits = 2;

    @Test
    public void testPhi(){
        RealMatrix alpha = new Array2DRowRealMatrix(new double[][] {
                {1.0, 3.0},
                {2.0, 4.0}
        });
        RealMatrix phiRM = OUPruneUtils.getPhiRM(tree.getNode(0), alpha);
        Double[] phiMat = new Double[] {2.3304025231210455, -1.5989838625030328, -1.0659892416686882, 0.7314186606180133};
        Double[] phiRMDouble = new Double[] {phiRM.getEntry(0,0), phiRM.getEntry(0,1), phiRM.getEntry(1,0), phiRM.getEntry(1,1)};
        Assert.assertArrayEquals(phiMat, phiRMDouble);
    }

    @Test
    public void testOmega(){
        RealVector theta = new ArrayRealVector(new double[] {0.5, 0.5});
        RealMatrix identity = MatrixUtils.createRealIdentityMatrix(nTraits);
        RealMatrix phiRM = new Array2DRowRealMatrix(new double[][]{
                {2.33040252312105, -1.59898386250303},
                {-1.06598924166869, 0.731418660618013}
        });
        RealVector omegaVec = OUPruneUtils.getOmegaVec(theta, phiRM, identity);
        Double[] omega = new Double[] {0.13429066969099002, 0.6672852905253386};
        Double[] omegaVecDouble = new Double[] {omegaVec.getEntry(0), omegaVec.getEntry(1)};
        Assert.assertArrayEquals(omega, omegaVecDouble);
    }

    @Test
    public void testAMatForOU(){
        RealMatrix inverseVarianceRM = new Array2DRowRealMatrix(new double[][]{
                {0.201407810205287, 0.215606098670895},
                {0.215606098670895, 0.64139862656716}
        });
        RealMatrix aRM = OUPruneUtils.getAMatForOU(inverseVarianceRM);
        Double[] aMat = new Double[] {-0.1007039051026435, -0.1078030493354475, -0.32069931328358};
        Double[] aRMDouble = new Double[] {aRM.getEntry(0,0), aRM.getEntry(0,1), aRM.getEntry(1,1)};
        Assert.assertArrayEquals(aMat, aRMDouble);
    }

    @Test
    public void testbVecForOU(){
        RealVector omega = new ArrayRealVector(new double[] {0.134290669690993, 0.667285290525338});
        RealMatrix invVariance = new Array2DRowRealMatrix(new double[][] {
                {0.201407810205287, 0.215606098670895,},
                {0.215606098670895, 0.64139862656716}
        });
        RealVector bVec = OUPruneUtils.getBVecForOU(invVariance, omega);
        Double[] b = new Double[] {0.17091796790410727, 0.4569497562513969};
        Double[] bVecDouble = new Double[] {bVec.getEntry(0), bVec.getEntry(1)};
        Assert.assertArrayEquals(b, bVecDouble);
    }

    @Test
    public void testCMatForOU(){

    }

    @Test
    public void testdVecForOU(){

    }

    @Test
    public void testEMatForOU(){
        RealMatrix phiMat = new Array2DRowRealMatrix(new double[][]{
                {11.7731179787193, -8.07801495941958},
                {-5.38534330627972, 3.69510301929973}});
        RealMatrix inverseMat = new Array2DRowRealMatrix(new double[][]{
                {0.128000066004857, 0.271161285503333,},
                {0.271161285503333, 0.599354305254683}
        });

        RealMatrix eRM = OUPruneUtils.getEMatForOU(phiMat,inverseMat);
        Double[] eMat = new Double[] {0.04666326455146064, -0.03531489040131186, -0.03201756321336657, 0.024230982265503354};
        Double[] eRMDouble = new Double[] {eRM.getEntry(0,0), eRM.getEntry(0,1), eRM.getEntry(1,0), eRM.getEntry(1,1)};
        Assert.assertArrayEquals(eMat, eRMDouble);
    }

    @Test
    public void testfForOU(){
        RealVector omega = new ArrayRealVector(new double[] {0.134290669690993, 0.667285290525338});
        RealMatrix invVariance = new Array2DRowRealMatrix(new double[][] {
                {0.201407810205287, 0.215606098670895,},
                {0.215606098670895, 0.64139862656716}
        });
        double detVariance = 12.0923805058647;
        double f = OUPruneUtils.getFforOU(omega, invVariance, detVariance, nTraits);
        Assert.assertEquals(-3.24809910801962,f, EPSILON);
    }

    @Test
    public void testPopulateUpperSigmaMatrix(){
        int nTraits = 3;
        RealMatrix sigmaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        double[] varValues = new double[]{3.6, 1.6, 2.5};
        double[] covValues = new double[] {0.2, 0.1, -0.5};
        OUPruneUtils.populateUpperSigmaMatrix(sigmaRM, varValues, covValues, nTraits);

        Double[] sigmaMat = new Double[] {3.6, 0.2, 0.1,
                                          0.0, 1.6, -0.5,
                                          0.0, 0.0, 2.5};
        Double[] sigmaRMDouble = new Double[] {sigmaRM.getEntry(0,0), sigmaRM.getEntry(0,1), sigmaRM.getEntry(0,2),
                                               sigmaRM.getEntry(1,0), sigmaRM.getEntry(1,1), sigmaRM.getEntry(1,2),
                                               sigmaRM.getEntry(2,0), sigmaRM.getEntry(2,1), sigmaRM.getEntry(2,2)};
        Assert.assertArrayEquals(sigmaMat, sigmaRMDouble);
    }

    @Test
    public void testPopulateSigmaMatrix () {
        int nTraits = 3;
        RealMatrix sigmaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        double[] sigmaValues = new double[]{3.6, 1.6, 2.5};
        double[] rhoValues = new double[] {0.2, 0.1, -0.5};
        OUPruneUtils.populateSigmaMatrix(sigmaRM, sigmaValues, rhoValues, nTraits);

        Double[] sigmaMat = new Double[] {3.6,                 0.48, 0.30000000000000004,
                                          0.48,                1.6,  -1.0,
                                          0.30000000000000004, -1.0,  2.5};
        Double[] sigmaRMDouble = new Double[] {sigmaRM.getEntry(0,0), sigmaRM.getEntry(0,1), sigmaRM.getEntry(0,2),
                                               sigmaRM.getEntry(1,0), sigmaRM.getEntry(1,1), sigmaRM.getEntry(1,2),
                                               sigmaRM.getEntry(2,0), sigmaRM.getEntry(2,1), sigmaRM.getEntry(2,2)};
        Assert.assertArrayEquals(sigmaMat, sigmaRMDouble);
    }

}
