package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.GeneralNodeMath;
import contraband.prunelikelihood.BMLikelihood;
import contraband.prunelikelihood.BMPruneShrinkageLikelihood;
import contraband.prunelikelihood.MorphologicalData;
import contraband.prunelikelihood.SigmaMatrix;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class BMLikelihoodTest {

    private final static double EPSILON = 1e-6;
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> contTraitData;
    private List<Double> popTraitData;
    private RealParameter correlation;
    private RealParameter rootValues;
    private RealParameter sigmasq;
    private RealParameter sigmaesq;
    private RealParameter delta;
    private RealParameter lambda;
    private final KeyRealParameter population = new KeyRealParameter();
    private final MorphologicalData morphData = new MorphologicalData();
    private final GeneralNodeMath nodeMath = new GeneralNodeMath();
    private final SigmaMatrix sigmaMatrix = new SigmaMatrix();
    private final SigmaMatrix sigmaEMatrix = new SigmaMatrix();
    private final KeyRealParameter contTrait = new KeyRealParameter();
    private final IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
    private final RealParameter colorValues = new RealParameter(new Double[] {1.0});
    private final RateCategoryClockModel lsc = new RateCategoryClockModel();

    @Test
    public void testBMLikelihood() {
        // tree
        treeStr = "(((t3:1.209461463,t4:1.209461463):0.6659547705,(t9:0.841016425,t10:0.841016425):1.034399809):1.561956365,(((t2:1.602817551,(t5:1.164343725,t6:1.164343725):0.4384738261):0.6643605462,t1:2.267178098):0.7120187616,(t7:1.115655119,t8:1.115655119):1.86354174):0.458175739);";
        spNames = "t3 t4 t9 t10 t2 t5 t6 t1 t7 t8";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 2;
        contTraitData =  Arrays.asList(
                0.326278727608277, -3.22668212260941,
                1.8164550628074, -1.71183724870188,
                -0.370085503473201, 1.81925405275285,
                0.665116986641999, -0.428821390843389,
                1.17377224776421, 4.22298205455098,
                3.59440970719762, 1.51483058860744,
                3.38137444987329, 3.63674837097173,
                -0.187743059073837, 3.68456953445085,
                -1.64759474375234, -0.743303344769609,
                -2.19534387982435, 1.10602125889508
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // BM model parameters
        sigmasq = new RealParameter(new Double[]{0.93818618860063, 1.74928729364035});
        correlation = new RealParameter(new Double[]{-0.117852694436286});
        sigmaMatrix.initByName("sigmasq", sigmasq, "correlation", correlation, "trait", morphData);

        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "tree", tree);

        // prune likelihood
        BMLikelihood pcm = new BMLikelihood();
        pcm.initByName("nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = pcm.calculateLogP();
        Assert.assertEquals(-39.588461762377, logP, EPSILON);
    }

    @Test
    public void testBMLikelihoodMatrix() {
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

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // BM model parameters
        rootValues = new RealParameter(new Double[]{2.0, -2.0, 1.5});
        sigmasq = new RealParameter(new Double[]{1.6, 2.5, 1.0});
        sigmaesq = new RealParameter(new Double[]{0.4, 0.9, 0.36});
        correlation = new RealParameter(new Double[]{0.5, -0.2, 0.1});
        sigmaMatrix.initByName("sigmasq", sigmasq, "correlation", correlation, "trait", morphData);
        sigmaEMatrix.initByName("sigmasq", sigmaesq, "correlation", correlation, "trait", morphData);
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "popMatrix", sigmaEMatrix, "tree", tree, "shareCorrelation", true, "rootValues", rootValues);

        Double[] sigmaValues = new Double[nTraits * nTraits];
        Double[] sigmaEValues = new Double[nTraits * nTraits];
        for (int i = 0; i < sigmaEValues.length; i ++) {
            sigmaValues[i] = sigmaMatrix.getSigmaMatrix()[i];
            sigmaEValues[i] = sigmaEMatrix.getSigmaMatrix()[i];
        }
        Assert.assertArrayEquals(new Double[]{1.6, 1.0, -0.2529822128134704, 1.0, 2.5, 0.158113883008419, -0.2529822128134704, 0.158113883008419, 1.0}, sigmaValues);
        Assert.assertArrayEquals(new Double[]{0.4, 0.3, -0.0758946638440411, 0.3, 0.9, 0.05692099788303082, -0.0758946638440411, 0.05692099788303082, 0.36}, sigmaEValues);


        // prune likelihood
        BMLikelihood pcm = new BMLikelihood();
        pcm.initByName("nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = pcm.calculateLogP();
        Assert.assertEquals(-18.2556630238236, logP, EPSILON);
    }

    @Test
    public void testBMLikelihoodStrictVariance() {
        treeStr = "((t1:0.506469612,t2:0.506469612):0.186771392,(t3:0.3878183895,t4:0.3878183895):0.3054226144);";
        spNames = "t1 t2 t3 t4";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 3;
        contTraitData =  Arrays.asList(
                -0.289918575728678, 1.34532611974715, -2.96529960582725,
                -0.712676267192549, 4.21710554156492, -1.28217844494824,
                0.477776523772291, 3.42975846405347, -2.20213758750652,
                2.33571507387484, 4.18494623610507, -3.99013430052531
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // BM model parameters
        rootValues = new RealParameter(new Double[]{0.0, 3.0, -3.0});
        sigmasq = new RealParameter(new Double[]{2.5});
        sigmaesq = new RealParameter(new Double[]{0.4});
        correlation = new RealParameter(new Double[]{0.5, -0.2, 0.1});
        sigmaMatrix.initByName("sigmasq", sigmasq, "correlation", correlation, "trait", morphData, "oneRateOnly", true);
        sigmaEMatrix.initByName("sigmasq", sigmaesq, "correlation", correlation, "trait", morphData, "oneRateOnly", true);
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "popMatrix", sigmaEMatrix, "tree", tree, "shareCorrelation", true, "rootValues", rootValues);

        // prune likelihood
        BMLikelihood pcm = new BMLikelihood();
        pcm.initByName("nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = pcm.calculateLogP();
        Assert.assertEquals(-18.6104281258402, logP, EPSILON);
    }

    @Test
    public void testBMLikelihoodWithShrinkage(){
        // tree
        treeStr = "((t4_2:4.1813368,t4_1:4.1813368):51.5649168,((t1_1:16.4777971,(t2_1:2.4163095,t3_1:2.4163095):14.0614876):7.5934926,t6_1:24.0712897):31.6749639):0.0;";
        spNames = "t4_2 t4_1 t1_1 t2_1 t3_1 t6_1";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // continuous trait values
        nTraits = 5;
        contTraitData = Arrays.asList(
                -1.49423749961236, -2.20443210635097, 13.0073187571802, -1.73557048290321, 0.362214907528917,
                1.61438393261534, -0.600823508278549, 6.76701308252399, 1.35749484021919, 0.0284817981254237,
                5.69171615806967, -1.00502666401073, 0.305337545316146, 0.330407167242628, 0.305460656456528,
                -1.76564797547353, 2.73908099097237, -10.2246819199159, 7.24123831044803, -0.143939688662232,
                -5.88260110808259, -0.304175058530707, -5.9742892461524, 0.674662257070058, -5.64650554205318,
                -11.0382115821499, 7.74032190916081, 2.29305900150846, -7.36164535234136, 7.52792561895395);
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);

        // population samples
        popTraitData = Arrays.asList(
                5.69171615806967, -1.00502666401073, 0.305337545316146, 0.330407167242628, 0.305460656456528,
                -1.76564797547353, 2.73908099097237, -10.2246819199159, 7.24123831044803, -0.143939688662232,
                -5.88260110808259, -0.304175058530707, -5.9742892461524, 0.674662257070058, -5.64650554205318,
                1.61438393261534, -0.600823508278549, 6.76701308252399, 1.35749484021919, 0.0284817981254237,
                -1.49423749961236, -2.20443210635097, 13.0073187571802, -1.73557048290321, 0.362214907528917,
                -11.0382115821499, 7.74032190916081, 2.29305900150846, -7.36164535234136, 7.52792561895395
        );
        String pop = "sp1 sp2 sp3 sp4 sp5 sp6";
        population.initByName("value", popTraitData, "keys", pop, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // shared evolutionary rate
        // correlations are estimated by shrinkage method
        sigmasq = new RealParameter(new Double[]{0.0467689});
        delta = new RealParameter(new Double[] {0.879396295242854});
        sigmaMatrix.initByName("sigmasq", sigmasq, "trait", morphData, "shrinkage", true, "delta", delta, "population", population, "oneRateOnly", true);

        // node math
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "tree", tree);

        // prune likelihood
        BMLikelihood pcm = new BMLikelihood();
        // here we subtract the likelihood at the root so that it matches mcmcTree
        pcm.initByName("nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc, "restrict", true);
        double logP = pcm.calculateLogP();
        Assert.assertEquals(-538.9563236764393, logP, EPSILON);
    }

    @Test
    public void testBMLikelihoodWithStandardizedData() {
        // tree
        treeStr = "((t2:0.8634660628,(t5:0.08432407594,t6:0.08432407594):1.623681925):1.067074613,((t3:0.2860224457,t4:0.2860224457):1.490933336,t1:1.776955782):0.9981248329);";
        spNames = "t2 t5 t6 t3 t4 t1";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // continuous trait values
        nTraits = 5;
        contTraitData = Arrays.asList(
                -0.79046984445445, -1.24550433961309, 1.00719121269088, -5.55377524207666, -0.190324276352273,
                -3.1456158063213, -4.40064210544044, -1.23674454689067, -4.92495067785601, 1.78818291334198,
                -3.1863298538836, -4.78828705170277, -1.89361406474235, -5.38105779594416, 3.01155546394615,
                -1.89832596392192, -1.4447079523426, 0.788657863129313, -4.36682104865865, -0.292610969178535,
                -4.41476753811221, -0.503260556700194, -1.07759730767263, -5.33637141004697, -1.92411328976331,
                1.24483751383358, -1.35069964377611, 0.0371698384792132, -6.31422173911139, -1.78833990891775
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);

        // population samples
        popTraitData = Arrays.asList(
                -1.49423749961236, -2.20443210635097, 13.0073187571802, -1.73557048290321, 0.362214907528917,
                1.61438393261534, -0.600823508278549, 6.76701308252399, 1.35749484021919, 0.0284817981254237,
                5.69171615806967, -1.00502666401073, 0.305337545316146, 0.330407167242628, 0.305460656456528,
                -1.76564797547353, 2.73908099097237, -10.2246819199159, 7.24123831044803, -0.143939688662232,
                -5.88260110808259, -0.304175058530707, -5.9742892461524, 0.674662257070058, -5.64650554205318,
                -11.0382115821499, 7.74032190916081, 2.29305900150846, -7.36164535234136, 7.52792561895395
        );
        String pop = "sp1 sp2 sp3 sp4 sp5 sp6";
        population.initByName("value", popTraitData, "keys", pop, "minordimension", nTraits);

        lambda = new RealParameter(new Double[]{0.604237758992233});
        morphData.initByName("traits", contTrait, "tree", tree, "population", population, "lambda", lambda);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // shared evolutionary rate
        // correlations are estimated by shrinkage method
        sigmasq = new RealParameter(new Double[]{2.5});
        delta = new RealParameter(new Double[] {0.879396295242854});
        sigmaMatrix.initByName("sigmasq", sigmasq, "trait", morphData, "shrinkage", true, "delta", delta, "population", population, "oneRateOnly", true);

        // accounting for intraspecific variance
        // after standardizing, all traits have the shared unit variance
        // the correlations are shared between species and individual samples
        sigmaesq = new RealParameter(new Double[]{1.0});
        sigmaEMatrix.initByName("sigmasq", sigmaesq, "correlation", correlation, "trait", morphData, "oneRateOnly", true);

        // node math
        rootValues = new RealParameter(new Double[]{-2.41413149877084, 0.55485848422132, 2.16888235336611, -4.6913954052587, 0.8582493776221});
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "popMatrix", sigmaEMatrix, "tree", tree, "rootValues", rootValues);

        // prune likelihood
        BMLikelihood pcm = new BMLikelihood();
        pcm.initByName("nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = pcm.calculateLogP();
        Assert.assertEquals(-56.1142534456042, logP, EPSILON);
    }


}
