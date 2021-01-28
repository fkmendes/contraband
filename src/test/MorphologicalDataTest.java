package test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.math.MatrixUtilsContra;
import contraband.prunelikelihood.MorphologicalData;
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

}
