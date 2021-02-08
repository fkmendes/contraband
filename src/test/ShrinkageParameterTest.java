package test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.prunelikelihood.MorphologicalData;
import contraband.prunelikelihood.ShrinkageParameter;
import contraband.utils.ShrinkageUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

public class ShrinkageParameterTest {
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> contTraitData;
    private final KeyRealParameter contTrait = new KeyRealParameter();
    private final MorphologicalData morphData = new MorphologicalData();
    private final ShrinkageParameter shrinkageParameter = new ShrinkageParameter();
    private int nSpecies;
    private double[] weight;
    private RealMatrix traitMat;
    private double delta;
    final static double EPSILON = 1e-7;

    @Test
    public void testShrinkageUtilsSixSpeciesFourTraits() {
        // initialize default weight vector
        nSpecies = 6;
        weight = new double [nSpecies];
        for (int i = 0; i < nSpecies; i++) {
            weight[i] = 1.0 / nSpecies;
        }

        // trait values
        traitMat = new Array2DRowRealMatrix(new double[][]{
                {1.43992069991891, -0.521710606400385, 1.10561637581434, -0.157612819068943},
                {0.348877438059058, 0.632563764185396, -0.157874694270553, -0.199347057660479},
                {-1.07819577781223, -0.466780487232386, -0.121903425587596, -0.735661003143472},
                { 1.08931843860522, -0.447764358092691, -0.200742989973232, -0.00327866708025785},
                {-1.14168814186063, -0.63123133899665, 0.01553830524637, -0.501053984275845},
                {-0.292510218946999, -1.56323463393269, -1.06283258902876, -0.75184508749932}
        });

        delta = ShrinkageUtils.estimateDelta(traitMat, weight);

        Assert.assertEquals(0.435026173002864, delta, EPSILON);
    }

    @Test
    public void testShrinkageUtilsFourSpeciesSixTraits(){
        // initialize default weight vector
        nSpecies = 4;
        weight = new double [nSpecies];
        for (int i = 0; i < nSpecies; i++) {
            weight[i] = 1.0 / nSpecies;
        }

        // trait values
        traitMat = new Array2DRowRealMatrix(new double[][]{
                {-0.420579569207178, 1.15032920533582, 0.0968521515822504, 1.8354220617026, 0.286678735881861, 0.836318736314219},
                {1.05096067831501, 0.158880815001583, 0.837912137759615, 1.47742100778513, -0.585856579077398, -0.366206661570247},
                {0.855288472381937, 0.524615741705375, 0.338739799123716, 0.0338961896720944, 1.39591803575264, -1.98905413542492},
                {0.745755760649241, -1.45262917873173, -2.46984489833647, 1.82339896004897, -0.900298166836283, 0.856360698349844}
        });

        delta = ShrinkageUtils.estimateDelta(traitMat, weight);

        Assert.assertEquals(0.676989791098336, delta, EPSILON);
    }

    @Test
    public void testShrinkageParameterSixSpeciesFourTraits(){
        // tree
        treeStr = "(t1:3.682105503,(t2:3.22892486,((t3:1.303684038,t4:1.303684038):1.056059251,(t5:1.773439676,t6:1.773439676):0.5863036128):0.8691815709):0.4531806428);";
        spNames = "t1 t2 t3 t4 t5 t6";
        tree = new TreeParser(treeStr, false, false, true, 0);
        nSpecies = 6;

        // trait values
        nTraits = 4;
        contTraitData =  Arrays.asList(
                1.43992069991891, -0.521710606400385, 1.10561637581434, -0.157612819068943,
                0.348877438059058, 0.632563764185396, -0.157874694270553, -0.199347057660479,
                -1.07819577781223, -0.466780487232386, -0.121903425587596, -0.735661003143472,
                1.08931843860522, -0.447764358092691, -0.200742989973232, -0.00327866708025785,
                -1.14168814186063, -0.63123133899665, 0.01553830524637, -0.501053984275845,
                -0.292510218946999, -1.56323463393269, -1.06283258902876, -0.75184508749932
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        Double[] weight = new Double [nSpecies];
        for (int i = 0; i < nSpecies; i++) {
            weight[i] = 1.0 / nSpecies;
        }
        RealParameter traitWeight = new RealParameter(weight);
        shrinkageParameter.initByName("trait", morphData, "weight", traitWeight);

        delta = shrinkageParameter.getDelta();
        Assert.assertEquals(0.435026173002864, delta, EPSILON);
    }

    @Test
    public void testShrinkageParameterFourSpeciesSixTraits(){
        // tree
        treeStr = "(((t1:0.2724306467,t2:0.2724306467):0.1437224097,t3:0.4161530564):0.745661126,t4:1.161814182);";
        spNames = "t1 t2 t3 t4";
        tree = new TreeParser(treeStr, false, false, true, 0);
        nSpecies = 4;

        // trait values
        nTraits = 6;
        contTraitData =  Arrays.asList(
                -0.420579569207178, 1.15032920533582, 0.0968521515822504, 1.8354220617026, 0.286678735881861, 0.836318736314219,
                1.05096067831501, 0.158880815001583, 0.837912137759615, 1.47742100778513, -0.585856579077398, -0.366206661570247,
                0.855288472381937, 0.524615741705375, 0.338739799123716, 0.0338961896720944, 1.39591803575264, -1.98905413542492,
                0.745755760649241, -1.45262917873173, -2.46984489833647, 1.82339896004897, -0.900298166836283, 0.856360698349844
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        Double[] weight = new Double [nSpecies];
        for (int i = 0; i < nSpecies; i++) {
            weight[i] = 1.0 / nSpecies;
        }
        RealParameter traitWeight = new RealParameter(weight);
        shrinkageParameter.initByName("trait", morphData, "weight", traitWeight);

        delta = shrinkageParameter.getDelta();
        Assert.assertEquals(0.676989791098336, delta, EPSILON);
    }
}
