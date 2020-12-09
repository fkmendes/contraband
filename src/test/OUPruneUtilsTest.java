package test;

import beast.evolution.tree.Tree;
import beast.util.TreeParser;
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

    /*
     * (1)
     */
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

    /*
     * (2)
     */
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

    /*
     * (3)
     */
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

    /*
     * (4)
     */
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

    /*
     * (5)
     */
    @Test
    public void testCMatForOU(){
        RealMatrix eMat = new Array2DRowRealMatrix(new double[][] {{0.239527487457355, -0.181275039198727}, {-0.164349514389387, 0.124380251933988}});
        RealMatrix phiMat = new Array2DRowRealMatrix(new double[][]{{2.33040252312105, -1.59898386250303}, {-1.06598924166869, 0.731418660618013}});

        RealMatrix cRM = OUPruneUtils.getCMatForOU(phiMat, eMat);

        Double[] cMat = new Double[] {-0.3757163513481894, 0.2577942667222091, -0.17688312929786393};
        Double[] cRMDouble = new Double[] {cRM.getEntry(0,0), cRM.getEntry(0,1), cRM.getEntry(1,1)};
        Assert.assertArrayEquals(cMat, cRMDouble);

    }

    /*
     * (6)
     */
    @Test
    public void testdVecForOU(){
        RealMatrix eMat = new Array2DRowRealMatrix(new double[][]{{0.239527487457355, -0.181275039198727}, {-0.164349514389387, 0.124380251933988}});
        RealVector omega = new ArrayRealVector(new double[] {0.134290669690993, 0.667285290525338});

        RealVector dVec = OUPruneUtils.getDVecForOU(eMat, omega);

        Double[] d = new Double[] {0.08879586049666546, -0.06092650619664564};
        Double[] dVecDouble = new Double[] {dVec.getEntry(0), dVec.getEntry(1)};
        Assert.assertArrayEquals(d, dVecDouble);
    }

    /*
     * (7)
     */
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

    /*
     * (8)
     */
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

    /*
     * (9)
     */
    @Test
    public void testOULMatForTips(){
        RealMatrix cMat = new Array2DRowRealMatrix(new double[][] {{-0.375716351348188, 0.257794266722209}, {0.257794266722209, -0.176883129297864}});

        RealMatrix lRM = OUPruneUtils.getLMatForOULeaf(cMat);

        Double[] lMat = new Double[] {-0.375716351348188, 0.257794266722209, -0.176883129297864};
        Double[] lRMDouble = new Double[] {lRM.getEntry(0, 0), lRM.getEntry(0, 1), lRM.getEntry(1, 1)};

        Assert.assertArrayEquals(lMat, lRMDouble);
    }

    /*
     * (10)
     */
    @Test
    public void testOUMVecForTips(){
        RealVector dVec = new ArrayRealVector(new double[]{0.0887958604966654, -0.0609265061966458});
        RealMatrix eMat = new Array2DRowRealMatrix(new double[][]{{0.239527487457355, -0.181275039198727}, {-0.164349514389387, 0.124380251933988}});
        RealVector traitVec = new ArrayRealVector(new double[]{4.53371989048144, 4.39037816144705});

        RealVector mVec = OUPruneUtils.getMVecForOULeaf(eMat, traitVec, dVec);

        Double[] m = new Double[] {0.3788804213855705, -0.25996482676851224};
        Double[] mVecDouble = new Double[] {mVec.getEntry(0), mVec.getEntry(1)};

        Assert.assertArrayEquals(m, mVecDouble);
    }

    /*
     * (11)
     */
    @Test
    public void testOURForTips(){
        RealVector traitVec = new ArrayRealVector(new double[]{4.53371989048144, 4.39037816144705});
        RealMatrix aMat = new Array2DRowRealMatrix(new double[][]{{-0.100703905102644, -0.107803049335448}, {-0.107803049335447, -0.32069931328358}});
        RealVector bVec = new ArrayRealVector(new double[]{0.170917967904107, 0.456949756251397});
        double f = -3.24809910801962;

        double r = OUPruneUtils.getRForOULeaf(aMat, traitVec, bVec, f);
        Assert.assertEquals(-13.0101512466986, r, EPSILON);
    }

    /*
     * (12)
     */
    @Test
    public void testOULMatForIntNode(){
        RealMatrix EAplusL_1 = new Array2DRowRealMatrix(new double[][]{{-0.428759247576708, 0.235066408456966}, {0.294189020642262, -0.161288815115095}});
        RealMatrix cMat = new Array2DRowRealMatrix(new double[][]{{-0.375716351348188, 0.257794266722209}, {0.257794266722209, -0.176883129297864}});
        RealMatrix eMat = new Array2DRowRealMatrix(new double[][]{{0.239527487457355, -0.181275039198727}, {-0.164349514389387, 0.124380251933988}});

        RealMatrix lRM = OUPruneUtils.getLMatForOUIntNode(cMat, eMat, EAplusL_1);

        Double[] lMat = new Double[] {-0.33938852692231414, 0.2328682684136375, -0.15978038776301526};
        Double[] lRMDouble = new Double[] {lRM.getEntry(0, 0), lRM.getEntry(0, 1), lRM.getEntry(1, 1)};

        Assert.assertArrayEquals(lMat, lRMDouble);
    }

    /*
     * (13)
     */
    @Test
    public void testOUMVecForIntNode(){
        RealMatrix EAplusL_1 = new Array2DRowRealMatrix(new double[][]{{-0.428759247576708, 0.235066408456966}, {0.294189020642262, -0.161288815115095}});
        RealVector dVec = new ArrayRealVector(new double[]{0.0887958604966654, -0.0609265061966458});
        RealVector bVec = new ArrayRealVector(new double[]{0.170917967904107, 0.456949756251397});
        RealVector mChild = new ArrayRealVector(new double[] {0.378880421385571, -0.259964826768511});

        RealVector mVec = OUPruneUtils.getMVecForOUIntNode(EAplusL_1, bVec, mChild, dVec);

        Double[] m = new Double[] {0.1835091624051842, -0.12591309810864493};
        Double[] mVecDouble = new Double[] {mVec.getEntry(0), mVec.getEntry(1)};

        Assert.assertArrayEquals(m, mVecDouble);
    }

    /*
     * (14)
     */
    @Test
    public void testOURForIntNode(){
        double rChild = -13.0101512466986;
        double f = -3.24809910801962;
        double logVNode = -0.152866886077582;
        RealMatrix AplusL_1 = new Array2DRowRealMatrix(new double[][]{{2.31907227633982, -0.699060988027452}, {-0.699060988027452, -2.22044210983387}});
        RealVector bVec = new ArrayRealVector(new double[]{0.170917967904107, 0.456949756251397});
        RealVector mChild = new ArrayRealVector(new double[] {0.378880421385571, -0.259964826768511});

        double r = OUPruneUtils.getRForOUIntNode(bVec, mChild, AplusL_1, f, rChild, 2, logVNode);
        Assert.assertEquals(-14.4597962945821, r, EPSILON);
    }

    /*
     * (15)
     */
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

    /*
     * (16)
     */
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

    /*
     * (17)
     */
    @Test
    public void testOmegaForJOU() {
        RealVector theta = new ArrayRealVector(new double[] {0.4, 0.6});
        RealMatrix identity = MatrixUtils.createRealIdentityMatrix(nTraits);
        RealMatrix phiRM = new Array2DRowRealMatrix(new double[][]{
                {2.33040252312105, -1.59898386250303},
                {-1.06598924166869, 0.731418660618013}
        });
        RealVector mj = new ArrayRealVector(new double[] {0.2, 0.8});
        RealVector iMinusPhiTheta = OUPruneUtils.getOmegaVec(theta, phiRM, identity);

        RealVector omegaVec = OUPruneUtils.getOmegaVec(mj, phiRM, iMinusPhiTheta);

        Double[] omega = new Double[] {-0.3858772771248161, 0.9594815804573406};
        Double[] omegaVecDouble = new Double[] {omegaVec.getEntry(0), omegaVec.getEntry(1)};

        Assert.assertArrayEquals(omega, omegaVecDouble);
    }

}
