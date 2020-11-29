package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.prunelikelihood.OUNodeMath;
import contraband.prunelikelihood.OUPruneLikelihood;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

public class OUPruneLikelihoodTest {
    final static double EPSILON = 1e-7;

    Tree tree;
    String treeStr;
    String spNames;
    int nTraits;
    List<Double> oneTraitValues;
    KeyRealParameter oneTraitData;

    OUNodeMath nodeMath;
    RealParameter rootValues;
    RealParameter sigma;
    RealParameter sigmae;
    RealParameter alpha;
    RealParameter theta;

    RealParameter colorValues;
    IntegerParameter colorAssignments;
    RateCategoryClockModel rcc;

    OUPruneLikelihood pcm;

    /*
     * (1) Tree with 3 taxa, 2 continuous traits, 1 optima
     */
    @Test
    public void testOUPruneLkSmallTreeTwoTraitsOneOpt() {
        treeStr = "((A:3.0058179,B:3.0058179):4.350951,C:7.3567689);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "A B C";

        nTraits = 2;
        oneTraitValues = Arrays.asList(
                0.131394584822684, -0.19269144948362,
                -0.65643500381027,-0.373259036803319,
                -0.17399894787897, 1.27056761078824
        );
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {0.2, 1.3});
        sigma = new RealParameter(new Double[] {1.0, 0.6, 1.0});
        sigmae = new RealParameter(new Double[] {1.0, 0.3, 1.0});
        alpha = new RealParameter(new Double[] {2.0, 0.0, 0.0, 2.0});
        theta = new RealParameter(new Double[] {0.5, 0.5});

        nodeMath = new OUNodeMath();
        colorAssignments = new IntegerParameter();
        colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0 });
        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta, "sigma", sigma, "sigmae", sigmae, "root", rootValues,
                "optNr", 1, "optAssign", colorAssignments);

        pcm = new OUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(-7.58111239313721, logP, EPSILON);
    }

    /*
     * (2) Tree with 4 taxa, 3 continuous traits, 1 optimum
     */
    @Test
    public void testOUPruneLkSmallTreeThreeTraitsOneOpt() {
        treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "sp1 sp2 sp3 sp4";

        nTraits = 3;
        oneTraitValues = Arrays.asList(
                2.90115170364898, 2.14268145872343, 3.2824264362529,
                2.62770319717086, 1.97031352144559, 4.28940778474205,
                3.96814434248716, 3.19423412470143, 3.49323574837925,
                8.88546166300816, 9.77485338099449, -3.65830173123419
        );
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        rootValues = new RealParameter(new Double[] {2.2, 1.3, 0.5});
        sigma = new RealParameter(new Double[] {0.35, 0.06, 0.06, 0.2, 0.4, 0.8});
        alpha = new RealParameter(new Double[] {1.8, -0.9, 1.2, -0.9, 1.0, 1.4, 1.2, 1.4, 1.6});
        theta = new RealParameter(new Double[] {5.0, 4.8, 1.2});
        nodeMath = new OUNodeMath();
        colorAssignments = new IntegerParameter();
        colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0});

        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta, "sigma", sigma, "sigmae", sigmae, "root", rootValues,
                "optNr", 1, "optAssign", colorAssignments);
        pcm = new OUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(-260.827821834006, logP, EPSILON);
    }


}
