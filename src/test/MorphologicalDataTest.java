package test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.math.GeneralNodeMath;
import contraband.math.MatrixUtilsContra;
import contraband.prunelikelihood.MorphologicalData;
import contraband.prunelikelihood.SigmaMatrix;
import org.apache.commons.lang3.ArrayUtils;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class MorphologicalDataTest {
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> contTraitData;
    private List<Double> popTraitData;
    private MorphologicalData morphData = new MorphologicalData();
    private final KeyRealParameter contTrait = new KeyRealParameter();
    private final KeyRealParameter population = new KeyRealParameter();
    private GeneralNodeMath nodeMath = new GeneralNodeMath();
    private final SigmaMatrix sigmaMatrix = new SigmaMatrix();
    private RealParameter sigmasq;
    private RealParameter delta;

    @Test
    public void testStandardizingData() {
        // tree
        treeStr = "((t2:0.5835094799,t3:0.5835094799):0.6852146462,t1:1.268724126);";
        spNames = "t2 t3 t1";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 3;
        contTraitData =  Arrays.asList(
                -0.151114078987955, 2.2129248459973, 1.66699822032105,
                0.726497005190457, 2.82635204778643, 3.97047716130522,
                2.67558373091379, -2.71828186018823, -3.48383931293939
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);

        // population samples
        popTraitData =  Arrays.asList(
                -1.34014546705207, 0.452826126748381, 1.5669866373813,
                2.39243664493101, 3.16853467161678, 0.784751288967007,
                -0.435195907092584, 0.626731799868515, -2.71660615533869,
                3.09845186877928, -0.0800885691840647, 3.22374078786359,
                3.52670145786031, 3.64869925640879, -1.64942346062992,
                -3.14424674902811, -0.104543044107213, -3.17031423924136
        );
        String popSpecies = "sp1 sp2 sp3 sp4 sp5 sp6";
        population.initByName("value", popTraitData, "keys", popSpecies , "minordimension", nTraits);

        RealParameter lambda = new RealParameter(new Double[]{0.729273516970115});
        morphData.initByName("traits", contTrait, "tree", tree, "population", population, "lambda", lambda);

        double[] standardizedData = morphData.getMorphologicalData();
        double[] t1 = new double[nTraits];
        MatrixUtilsContra.getMatrixRow(standardizedData, 0, nTraits, t1);
        double[] t2 = new double[nTraits];
        MatrixUtilsContra.getMatrixRow(standardizedData, 1, nTraits, t2);
        double[] t3 = new double[nTraits];
        MatrixUtilsContra.getMatrixRow(standardizedData, 2, nTraits, t3);
        // species t1
        Assert.assertArrayEquals(new Double[] { 1.025746456929693, -1.1517236447753618, -1.356609049223971}, ArrayUtils.toObject(t1));
        // species t2
        Assert.assertArrayEquals(new Double[] { -0.05793305188813881, 0.9376061793200064, 0.6491300739182903 }, ArrayUtils.toObject(t2));
        // species t3
        Assert.assertArrayEquals(new Double[] { 0.2785193079304739, 1.1975124910959332, 1.5461061096469335 }, ArrayUtils.toObject(t3));
    }

    @Test
    public void testTransformData() {
        // tree
        treeStr = "(((t4:0.3029120562,t5:0.3029120562):0.1726641058,t3:0.475576162):1.976889551,(t1:1.296371029,t2:2.035157287):0.4173084261);";
        spNames = "t4 t5 t3 t1 t2";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 3;
        contTraitData =  Arrays.asList(
                -2.58346070817778, -1.55079811333479, 5.6255436771746,
                -2.31424616039559, -1.23120551561249, 6.24842157812095,
                -1.59490654203882, -2.74721796020631, 7.27376543652255,
                -2.09305102331369, 1.70943420509975, 0.24750026425326,
                -0.333792899540118, 0.701543030250452, 1.55368477226022
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);

        // population samples
        popTraitData =  Arrays.asList(
                10.6890174876092, -3.52322456753132, 1.00287532421597, 1.3608173901569, -2.0037573466527, 7.09534979105388, -3.37399440946941, -3.00582928516136, 4.8214048794019, 8.91203567871889, -2.90221314494706, -4.73280230027901, -0.624935005660081, -1.3861790153386, 2.01716602752887, -7.19366892895753, -10.6876285148211, 6.5793185171151, -7.53145456709179, 8.44178549772067, -4.50513153957943, -0.796723866385605, 7.34302049395171, -4.50026866975631
        );
        String popSpecies = "sp1 sp2 sp3 sp4 sp5 sp6 sp7 sp8";
        population.initByName("value", popTraitData, "keys", popSpecies , "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree, "transform", true);

        // BM model parameters
        delta = new RealParameter(new Double[]{0.676132776830971});
        sigmasq = new RealParameter(new Double[]{2.77639276177938, 2.82578533469758, 4.54723563720239});
        sigmaMatrix.initByName("sigmasq", sigmasq, "trait", morphData, "delta", delta, "shrinkage", true, "population", population);

        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "tree", tree);

        morphData.transformTraitData(nodeMath.getTraitRateMatrix());
        double[] transformedData = morphData.getMorphologicalData();
        double[] t1 = new double[nTraits];
        MatrixUtilsContra.getMatrixRow(transformedData, 0, nTraits, t1);
        double[] t2 = new double[nTraits];
        MatrixUtilsContra.getMatrixRow(transformedData, 1, nTraits, t2);
        double[] t3 = new double[nTraits];
        MatrixUtilsContra.getMatrixRow(transformedData, 2, nTraits, t3);
        double[] t4 = new double[nTraits];
        MatrixUtilsContra.getMatrixRow(transformedData, 3, nTraits, t4);
        double[] t5 = new double[nTraits];
        MatrixUtilsContra.getMatrixRow(transformedData, 4, nTraits, t5);

        // species t1
        Assert.assertArrayEquals(new Double[] { -1.256143813303816, 0.9533069164718878, 0.260715297698597 }, ArrayUtils.toObject(t1));
        // species t2
        Assert.assertArrayEquals(new Double[] { -0.20032568772176615, 0.4075325993599785, 0.8356964501754756 }, ArrayUtils.toObject(t2));
        // species t3
        Assert.assertArrayEquals(new Double[] { -0.9571825833505195, -1.685954626381353, 3.041769729442216 }, ArrayUtils.toObject(t3));
        // species t4
        Assert.assertArrayEquals(new Double[] { -1.550463007993594, -1.0039576231191356, 2.3672975030234613 }, ArrayUtils.toObject(t4));
        // species t5
        Assert.assertArrayEquals(new Double[] { -1.3888939946818242, -0.8052281439696123, 2.728127406177297 }, ArrayUtils.toObject(t5));
    }

}
