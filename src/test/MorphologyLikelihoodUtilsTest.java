package test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.math.GeneralNodeMath;
import contraband.math.MatrixUtilsContra;
import contraband.prunelikelihood.MorphologicalData;
import contraband.prunelikelihood.SigmaMatrix;
import contraband.utils.MorphologyLikelihoodUtils;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class MorphologyLikelihoodUtilsTest {
    private final static double EPSILON = 1e-6;
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> contTraitData;
    MorphologicalData morphData = new MorphologicalData();
    private RealParameter rootValues;
    private RealParameter sigmasq;
    private RealParameter sigmaesq;
    private RealParameter correlation;
    private final GeneralNodeMath nodeMath = new GeneralNodeMath();
    private final SigmaMatrix sigmaMatrix = new SigmaMatrix();
    private final SigmaMatrix sigmaEMatrix = new SigmaMatrix();
    private final KeyRealParameter contTrait = new KeyRealParameter();

    @Test
    public void testMatrixVersion() {
        // tree
        treeStr = "((t3:0.276070382,t4:0.276070382):0.7116642276,(t1:0.4249275045,t2:0.4249275045):0.562807105);";
        spNames = "t3 t4 t1 t2";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 3;
        contTraitData =  Arrays.asList(
                1.306587759209, -3.83501774599227, 2.17725305921907,
                1.10688487066549, -0.962594500743626, 3.25013136469802,
                2.42388899371546, -1.49001389065793, 2.249184569853,
                4.06586793293621, -0.629432336665211, 0.887441181665878
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);
        double[] traitValuesArr = new double[]{2.42388899371546, -1.49001389065793, 2.249184569853,4.06586793293621, -0.629432336665211, 0.887441181665878, 1.306587759209, -3.83501774599227, 2.17725305921907, 1.10688487066549, -0.962594500743626, 3.25013136469802};

        // BM model parameters
        rootValues = new RealParameter(new Double[]{2.0, -2.0, 1.5});
        sigmasq = new RealParameter(new Double[]{1.6, 2.5, 1.0});
        sigmaesq = new RealParameter(new Double[]{0.4, 0.9, 0.36});
        correlation = new RealParameter(new Double[]{0.5, -0.2, 0.1});
        sigmaMatrix.initByName("sigmasq", sigmasq, "correlation", correlation, "trait", morphData);
        sigmaEMatrix.initByName("sigmasq", sigmaesq, "correlation", correlation, "trait", morphData);
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "popMatrix", sigmaEMatrix, "tree", tree, "shareCorrelation", true, "rootValues", rootValues);

        // populate trait rate matrix
        nodeMath.updateSigmaMatrix();
        nodeMath.operateOnTraitRateMatrix();
        nodeMath.operateOnInvTraitRateMatrix();

        // iterate each node and calculate A C E f
        for(int nodeIdx = 0; nodeIdx < tree.getNodeCount(); nodeIdx ++) {
            double branchLength = tree.getNode(nodeIdx).getLength();
            if(nodeIdx < tree.getLeafNodeCount()) {
                MorphologyLikelihoodUtils.populateACEfMatrixForTips(nodeMath, branchLength, nTraits, nodeIdx);
            } else {
                MorphologyLikelihoodUtils.populateACEfMatrixForIntNodes(nodeMath, branchLength, nTraits, nodeIdx);
            }
        }

        // iterate each node and calculate L m r
        for(int nodeIdx = 0; nodeIdx < tree.getNodeCount(); nodeIdx ++) {
            if(!tree.getNode(nodeIdx).isRoot()) {
            if(nodeIdx < tree.getLeafNodeCount()) {
                MorphologyLikelihoodUtils.populateLmrMatrixForTip(nodeMath, traitValuesArr, nTraits, nodeIdx);
            } else {
                MorphologyLikelihoodUtils.populateLmrMatrixForIntNode(nodeMath,nTraits, nodeIdx);
            }

                double[] lMat = nodeMath.getLMatForNode(nodeIdx).clone();
                double[] mVec = nodeMath.getMVecForNode(nodeIdx).clone();
                double r = nodeMath.getRForNode(nodeIdx);

                int pdx = tree.getNode(nodeIdx).getParent().getNr();
                double[] lp = nodeMath.getLMatForNode(pdx);
                double[] mp = nodeMath.getMVecForNode(pdx);
                double rp = nodeMath.getRForNode(pdx);

                double[] resl = new double[nTraits * nTraits];
                double[] resm = new double[nTraits];
                MatrixUtilsContra.vectorAdd(lMat, lp, resl);
                MatrixUtilsContra.vectorAdd(mVec, mp, resm);
                nodeMath.setLForNode(pdx, resl);
                nodeMath.setMVecForNode(pdx, resm);
                nodeMath.setRForNode(pdx, rp + r);
            }
        }

        // test intermediate parameters, A C E f and L m r
        int node1 = 0;
        Double[] aMat1 = new Double[nTraits * nTraits];
        Double[] cMat1 = new Double[nTraits * nTraits];
        Double[] eMat1 = new Double[nTraits * nTraits];
        double f1 = nodeMath.getfForNode(node1);
        Double[] lMat1 = new Double[nTraits * nTraits];
        Double[] mVec1 = new Double[nTraits];
        double r1 = nodeMath.getRForNode(node1);

        int node2 = 5;
        Double[] aMat2 = new Double[nTraits * nTraits];
        Double[] cMat2 = new Double[nTraits * nTraits];
        Double[] eMat2 = new Double[nTraits * nTraits];
        double f2 = nodeMath.getfForNode(node2);
        Double[] lMat2 = new Double[nTraits * nTraits];
        Double[] mVec2 = new Double[nTraits];
        double r2 = nodeMath.getRForNode(node2);

        for (int i = 0; i < nTraits * nTraits; i ++) {
            aMat1[i] = nodeMath.getAMatForNode(node1)[i];
            aMat2[i] = nodeMath.getAMatForNode(node2)[i];

            cMat1[i] = nodeMath.getCMatForNode(node1)[i];
            cMat2[i] = nodeMath.getCMatForNode(node2)[i];

            eMat1[i] = nodeMath.getEMatForNode(node1)[i];
            eMat2[i] = nodeMath.getEMatForNode(node2)[i];

            lMat1[i] = nodeMath.getLMatForNode(node1)[i];
            lMat2[i] = nodeMath.getLMatForNode(node2)[i];
        }

        for (int i = 0; i < nTraits; i ++) {
            mVec1[i] = nodeMath.getMVecForNode(node1)[i];
            mVec2[i] = nodeMath.getMVecForNode(node2)[i];
        }

        Assert.assertArrayEquals(new Double[]{-0.6716360815648698, 0.26064967365653974, -0.19813621157073422, 0.26064967365653974, -0.3585277135662245, 0.11758749522332192, -0.19813621157073422, 0.1175874952233219, -0.7018870963188379}, aMat1);
        Assert.assertArrayEquals(new Double[]{-0.6716360815648698, 0.26064967365653974, -0.19813621157073422, 0.26064967365653974, -0.3585277135662245, 0.11758749522332192, -0.19813621157073422, 0.1175874952233219, -0.7018870963188379}, cMat1);
        Assert.assertArrayEquals(new Double[]{1.3432721631297395, -0.5212993473130795, 0.39627242314146843, -0.5212993473130795, 0.717055427132449, -0.23517499044664383, 0.39627242314146843, -0.2351749904466438, 1.4037741926376759}, eMat1);
        Assert.assertArrayEquals(new Double[]{-0.6716360815648698, 0.26064967365653974, -0.19813621157073422, 0.26064967365653974, -0.3585277135662245, 0.11758749522332192, -0.19813621157073422, 0.1175874952233219, -0.7018870963188379}, lMat1);
        Assert.assertArrayEquals(new Double[]{4.923975700049984, -2.8609462569101485, 4.46828162110492}, mVec1);

        Assert.assertArrayEquals(new Double[]{-0.8083821814831524, 0.3396838257545367, -0.2582150417557933, 0.33968382575453676, -0.5016868811143927, 0.1652576267237077, -0.2582150417557933, 0.1652576267237077, -0.9798571896765483}, aMat2);
        Assert.assertArrayEquals(new Double[]{-0.8083821814831524, 0.3396838257545367, -0.2582150417557933, 0.33968382575453676, -0.5016868811143927, 0.1652576267237077, -0.2582150417557933, 0.1652576267237077, -0.9798571896765483}, cMat2);
        Assert.assertArrayEquals(new Double[]{1.6167643629663049, -0.6793676515090734, 0.5164300835115866, -0.6793676515090735, 1.0033737622287855, -0.3305152534474154, 0.5164300835115866, -0.3305152534474154, 1.9597143793530967}, eMat2);
        Assert.assertArrayEquals(new Double[]{-0.5043263094152102, 0.2052636823310822, -0.15603383583641356, 0.20526368233108225, -0.29502298731210885, 0.09694450521473906, -0.15603383583641356, 0.09694450521473909, -0.5769741539210209}, lMat2);
        Assert.assertArrayEquals(new Double[]{4.197420244535865, -2.2614753931552274, 3.027842321897972}, mVec2);

        Assert.assertEquals(-2.8202176061436, f1, EPSILON);
        Assert.assertEquals(-2.39490404560275, f2, EPSILON);

        Assert.assertEquals(-15.9442227280648, r1, EPSILON);
        Assert.assertEquals(-18.4178304735558, r2, EPSILON);
    }
}
