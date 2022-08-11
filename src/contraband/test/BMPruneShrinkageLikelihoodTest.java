package contraband.test;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeParser;
import contraband.prunelikelihood.BMPruneLikelihood;
import contraband.clock.RateCategoryClockModel;
import contraband.math.NodeMath;
import contraband.prunelikelihood.BMPruneShrinkageLikelihood;
import org.junit.Assert;
import org.junit.Test;
import java.util.Arrays;
import java.util.List;

/*
 * This class contains unit tests for BMPruneShrinkageLikelihood.class
 * (using shrinkage method)
 */

public class BMPruneShrinkageLikelihoodTest {
    final static double EPSILON = 1e-6;
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> data;
    private RealParameter sigmasq;
    private final NodeMath nodeMath = new NodeMath();
    private final RealParameter traitValues = new RealParameter();
    private final IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
    private RealParameter colorValues = new RealParameter(new Double[] {1.0});
    private final RateCategoryClockModel lsc = new RateCategoryClockModel();

    /*
     * tree with 6 species and each species has 5 continuous traits
     * t4_1 and t6_1 are fossils
     * the rest of the taxa are extant species
     */
    @Test
    public void testBMPruneShrinkageLikelihood6Species5Traits() {
        // tree
        treeStr = "((t4_2:9.1379169,t4_1:8.470416972):33.6333243,((t1_1:16.0044120,(t2_1:5.7809186,t3_1:5.7809186):10.2234934):3.4863988,t6_1:14.863633601):23.2804304):0.0;";
        spNames = "t4_2 t4_1 t1_1 t2_1 t3_1 t6_1";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 5;
        data = Arrays.asList(
                -1.49423749961236, -2.20443210635097, 13.0073187571802, -1.73557048290321, 0.362214907528917,
                -1.01313598636176, -1.45295055296761, 8.20086495770802, 0.459868047614318, 0.494318850519285,
                5.69171615806967, -1.00502666401073, 0.305337545316146, 0.330407167242628, 0.305460656456528,
                -1.76564797547353, 2.73908099097237, -10.2246819199159, 7.24123831044803, -0.143939688662232,
                -5.88260110808259, -0.304175058530707, -5.9742892461524, 0.674662257070058, -5.64650554205318,
                -3.25522309907749, 5.68185395502633, 1.38164443830931, -4.82948242309325, 4.23511851679411
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // node math
        sigmasq = new RealParameter(new Double[]{1.0479483});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "shrinkage", true, "oneRateOnly", true);

        // likelihood
        BMPruneShrinkageLikelihood PCM4 = new  BMPruneShrinkageLikelihood();
        PCM4.initByName("nodeMath", nodeMath, "tree", tree, "traits",  traitValues, "branchRateModel", lsc, "delta", 0.963523553937718);
        double lnLk = PCM4.calculateLogP();

        Assert.assertEquals(-76.65483454526472, lnLk, EPSILON);
    }

    /*
     * tree with 3 species and each species has 2 continuous traits
     * taxa A and B are extant species
     * taxa C is a fossil
     */
    @Test
    public void testBMPruneShrinkageLikelihood3Species2Traits() {
        // tree
        treeStr = "((A:12.4420263,B:12.4420263):42.9258211,C:43.5702874);";
        spNames = "A B C";
        tree = new TreeParser(treeStr, false, false, true, 0);

        //trait values
        nTraits = 2;
        data = Arrays.asList(
                1.0, 2.0,
                3.0, 5.0,
                2.0, 4.0
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // node math
        sigmasq = new RealParameter(new Double[]{0.1543038});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "shrinkage", true, "oneRateOnly", true);

        // likelihood
        BMPruneShrinkageLikelihood PCM2 = new BMPruneShrinkageLikelihood();
        PCM2.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc, "delta", 0.25925925925926);
        double lnLk = PCM2.calculateLogP();

        Assert.assertEquals(-8.128457538378111, lnLk, EPSILON);
    }

    /*
     * tree with 6 species and each species has 5 continuous traits
     * all the taxa are extant species
     */
    @Test
    public void testBMPruneShrinkageLikelihood6Species5TraitsUltrametricTree() {
        // tree
        treeStr = "((t4_2:4.1813368,t4_1:4.1813368):51.5649168,((t1_1:16.4777971,(t2_1:2.4163095,t3_1:2.4163095):14.0614876):7.5934926,t6_1:24.0712897):31.6749639):0.0;";
        spNames = "t4_2 t4_1 t1_1 t2_1 t3_1 t6_1";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 5;
        data = Arrays.asList(
                -1.49423749961236, -2.20443210635097, 13.0073187571802, -1.73557048290321, 0.362214907528917,
                1.61438393261534, -0.600823508278549, 6.76701308252399, 1.35749484021919, 0.0284817981254237,
                5.69171615806967, -1.00502666401073, 0.305337545316146, 0.330407167242628, 0.305460656456528,
                -1.76564797547353, 2.73908099097237, -10.2246819199159, 7.24123831044803, -0.143939688662232,
                -5.88260110808259, -0.304175058530707, -5.9742892461524, 0.674662257070058, -5.64650554205318,
                -11.0382115821499, 7.74032190916081, 2.29305900150846, -7.36164535234136, 7.52792561895395);
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // node math
        sigmasq = new RealParameter(new Double[]{0.0467689});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "shrinkage", true, "oneRateOnly", true);

        // likelihood
        BMPruneShrinkageLikelihood PCM3 = new BMPruneShrinkageLikelihood();
        PCM3.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc, "delta", 0.879396295242854);
        double lnLk = PCM3.calculateLogP();

        Assert.assertEquals(-538.9563236764393, lnLk, EPSILON);
    }

    /*
     * tree with 3 species and each species has 2 continuous traits
     * the original data are used to estimate the measurement error of the population
     */
    @Test
    public void testBMPruneShrinkageLikelihood3Species2TraitsPopSE() {
        // tree
        treeStr = "((A:23.3161955,B:23.3161955):16.0514844,C:39.3676799);";
        spNames = "A B C";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 2;
        // original data
        data = Arrays.asList(
                -2.62762948691895, -1.56292164859448,
                -1.50846427625826, -1.59482814741543,
                -0.226074849617958, -2.11000367246907
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // branch rate model
        colorValues = new RealParameter(new Double[]{0.3465491});
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // node math
        sigmasq = new RealParameter(new Double[]{1.0});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "shrinkage", true, "oneRateOnly", true);

        // likelihood
        RealParameter popVar = new RealParameter(new Double[]{0.3});
        BMPruneShrinkageLikelihood PCM6 = new BMPruneShrinkageLikelihood();
        PCM6.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc, "delta", 0.303252482035699, "includePopVar", true, "popVar", popVar);
        double lnLk = PCM6.calculateLogP();
        Assert.assertEquals(-9.73103796742531, lnLk, EPSILON);
    }

    /*
     * tree with 5 species and each species has 9 continuous traits
     *
     */
    @Test
    public void testBMPruneWithWithoutShrinkageLikelihood() {
        // number of traits
        nTraits = 9;

        // species names
        spNames = "t2 t4 t5 t3 t1";

        // original data
        data = Arrays.asList(
                0.996206925293026, -2.91270268270845, 1.35994187362101, 0.102975767496277, -2.43052504782334, 2.98602562445202, -1.03975429308448, -10.8901344038457, -1.30156469246422,
                -2.65052276370749, -1.5792798388905, 6.82985203851052, 1.79880717599135, -0.867723116858107, 7.49844240046301, 4.50781790071978, 0.520102801698764, -0.00971374138577277,
                -3.93588311090379, -6.4085861951559, 9.63880711538996, 1.31978890202522, -0.637507941646737, 8.33514902203189, 8.4164152934786, 0.429643092693121, -1.47520049052425,
                -5.59242339046684, -5.36037547409755, 8.7939796440525, 3.92494613908588, 0.0934316079593845, 8.15176984162606, 6.70696946959889, -2.11120755504519, 0.74040432381528,
                -2.53243812537148, 5.72854673050833, 5.53576146755777, 1.16992957583461, -0.550563880149175, 2.86099545245125, 11.752265679971, -1.08417023425239, -4.43344316554345);
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // (1) BMPruneLikelihood using shrinkage
        // res.tree1, this is a sampled tree from mcmcTree
        String treeStr1 = "((t2:21.5238014,((t4:2.7950731,t5:2.7950731):9.1838471,t3:11.9789202):9.5448812):19.9350992,t1:41.4589006);";
        TreeParser tree1 = new TreeParser(treeStr1, false, false, true, 0);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree1);

        // node math
        NodeMath nodeMath1 = new NodeMath();
        RealParameter sgimasq1 = new RealParameter(new Double[]{0.2125345});
        nodeMath1.initByName("traits", traitValues, "sigmasq", sgimasq1, "oneRateOnly", true, "shrinkage", true);

        // likelihood
        BMPruneShrinkageLikelihood PCM1 = new BMPruneShrinkageLikelihood();
        PCM1.initByName("nodeMath", nodeMath1, "tree", tree1, "traits", traitValues, "branchRateModel", lsc, "delta", 0.582249579199718);
        double lnLk1 = PCM1.calculateLogP();
        // expected: -101.547
        Assert.assertEquals(-101.54729551277926, lnLk1, EPSILON);

        // 2. PCM likelihood without using shrinkage
        // res.tree2, this is a sampled tree from mcmcTree
        String treeStr2 = "((t2:21.0771224,((t4:1.2011247,t5:1.2011247):12.1000181,t3:13.3011428):7.7759796):16.2753358,t1:37.3524582);";
        TreeParser tree2 = new TreeParser(treeStr2, false, false, true, 0);

        // rate matrix without shrinkage, i.e. true values when simulating data
        // Note: Sigma <- sigma * t(sigma)
        RealParameter covariance = new RealParameter(new Double[] {
                0.350677950781843, -0.236395703863967, -0.533911852993423, -0.108987457329148, -0.51300222415362, -0.36444561464557, -0.312519633924628, -0.843346655371207,
                0.418893522247262, 0.0766865589182576, -0.569068468505974, 0.626907460534963, 0.213232110738407, -0.147535741496136, 0.447562830522511,
                0.439066743766521, 0.410790540818744, 0.344320127050045, 0.276958826768794, -0.0309558552883608, -0.152981331869958,
                -0.190235500506191, -0.347353489460256, -0.103958639009799, -0.632698175784132, 1.08447798266926,
                0.603980999152495, -0.561554291601675, -0.201442417649538, -0.233327676811609,
                0.389982559168159, -0.0416845332359146, 0.126659256997377,
                -0.0142733776743515, -0.021435228645658,
                0.684301142007229
        });
        RealParameter variance = new RealParameter(new Double[] {
                2.37809895230671, 3.92444342735533, 3.02812661102414,
                3.06556360045831, 2.81005713192123, 2.30077876608194,
                4.64631271741818, 3.15615698487533, 1.50684016606037});

        // branch rate model
        RateCategoryClockModel lsc2 = new RateCategoryClockModel();
        IntegerParameter colorAssignments2 = new IntegerParameter(new Integer[]{0});
        RealParameter colorValues2 = new RealParameter(new Double[]{0.2061328});
        lsc2.initByName("nCat", 1, "rateCatAssign", colorAssignments2, "rates", colorValues2, "tree", tree2);

        // node math
        NodeMath nodeMath2 = new NodeMath();
        nodeMath2.initByName("traits", traitValues, "sigmasq", variance, "covariance", covariance, "upperMatrix", true);

        // likelihood
        BMPruneLikelihood PCM2 = new BMPruneLikelihood();
        PCM2.initByName("nodeMath", nodeMath2, "tree", tree2, "traits", traitValues, "branchRateModel", lsc2, "includeRoot", true);
        double lnLk2 = PCM2.calculateLogP();
        // expected: -102.655

        Assert.assertEquals(-102.65489905306194, lnLk2, EPSILON);
    }
}
