package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.NodeMath;
import contraband.prunelikelihood.BMPruneLikelihood;
import contraband.prunelikelihood.BMPruneShrinkageLikelihood;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class BMPruneLikelihoodWithMissingData {
    private final static double EPSILON = 1e-6;
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> data;
    private RealParameter correlation;
    private RealParameter sigmasq;
    private RealParameter rootValues;
    private final NodeMath nodeMath = new NodeMath();
    private final KeyRealParameter traitValues = new KeyRealParameter();
    private final IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
    private final RealParameter colorValues = new RealParameter(new Double[] {1.0});
    private final RateCategoryClockModel lsc = new RateCategoryClockModel();
/*
    @Test
    public void testBMPruneLikelihoodTwoOfFiveSpeciesMissingData() {
        // tree
        treeStr = "(((t4:4.157956267,(t5:3.645264322,t6:3.645264322):0.5126919446):3.853008658,(t2:5.632320957,t3:5.632320957):2.378643968):4.138629419,(((t8:2.935100956,(t9:1.420145405,t10:1.420145405):1.514955551):0.5094521248,t7:3.444553081):7.252722915,t1:10.697276):1.452318347);";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 3;
        data = Arrays.asList(
              -3.56841516408629, -5.93547758767035, -0.274921026230003,
              -3.1604638705923, 6.96623701978737, 2.44891129101194,
              14.3661798575386, 1.90468950651659, -17.2980486429996,
              12.1977107783458, -1.42445780496105, -16.8458053660924,
              6.93338326603582, -6.98429889078827, -12.7589216423731,
              0.410354347208799, -8.34121610664249, -4.56966693197476
        );
        spNames = "t6 t3 t9 t10 t7 t1";
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // BM model parameters
        sigmasq = new RealParameter(new Double[]{3.23149028594509,4.44002933545446,3.07371621624534});
        correlation = new RealParameter(new Double[]{0.407831756026716, -0.929925634986149, -0.0727691475140235});
        // the MLE has not been implemented yet. So we have to specify the root so that it matches mvMORPH.
        rootValues = new RealParameter(new Double[] {1.88992288820042, -2.81680659985333, -5.39528907654654});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation, "tree", tree, "rootValues", rootValues);

        // prune likelihood
        BMPruneLikelihood fpcm1 = new BMPruneLikelihood();
        fpcm1.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lik1 = fpcm1.calculateLogP();

        Assert.assertEquals(-43.7932520448502, lik1, EPSILON);
    }
*/

    @Test
    public void testBMPruneShrinkageLikelihoodTwoOfFiveSpeciesMissingData() {
        // tree
        treeStr = "((t4_2:4.1813368,t4_1:4.1813368):51.5649168,(((t1_1:10.0,t8_1:10.0):6.4777971,(t2_1:2.4163095,t3_1:2.4163095):14.0614876):7.5934926, (t6_1:20.0, t7_1:20.0):4.0712897):31.6749639):0.0;";
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
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "shrinkage", true, "oneRateOnly", true, "tree", tree);

        // likelihood
        BMPruneShrinkageLikelihood PCM3 = new BMPruneShrinkageLikelihood();
        PCM3.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc, "delta", 0.879396295242854);
        double lnLk = PCM3.calculateLogP();

        Assert.assertEquals(-538.9563236764393, lnLk, EPSILON);
    }

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

}
