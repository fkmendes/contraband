package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.prunelikelihood.OUNodeMath;
import contraband.prunelikelihood.DOUPruneLikelihood;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class DOUPruneLikelihoodTest {
    final static double EPSILON = 1e-7;

    Tree tree;
    String treeStr;
    String spNames;
    int nTraits;
    List<Double> oneTraitValues;
    KeyRealParameter oneTraitData;

    OUNodeMath nodeMath;
    RealParameter rootValues;
    RealParameter sigmasq;
    RealParameter covariance;
    RealParameter popVar;
    RealParameter popCov;
    RealParameter alpha;
    RealParameter dAlpha;
    RealParameter theta;
    IntegerParameter colorAssignments;
    DOUPruneLikelihood pcm;

    /*
     * (1) Ultrametric tree with 3 taxa, 2 continuous traits, 1 optima
     */
    @Test
    public void testDOUPruneLkSmallTreeTwoTraitsOneOpt() {
        treeStr = "((A:3.0,B:3.0):4.0,C:7.0);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "A B C";

        nTraits = 2;
        oneTraitValues = Arrays.asList(
                -8.35873276897194, -8.75604744861027,
                -10.4625464429157, -10.1448377936672,
                0.381165633761954, 2.85191670820294
        );
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {0.3, 0.4});
        sigmasq = new RealParameter(new Double[] {1.0, 1.0});
        covariance = new RealParameter(new Double[] {0.6});
        popVar = new RealParameter(new Double[] {1.0, 1.0});
        popCov = new RealParameter(new Double[] {0.4});
        alpha = new RealParameter(new Double[] {-0.2, 0.1, 0.1, -0.2});
        dAlpha = new RealParameter(new Double[] {0.1, -0.5, -0.5, 0.1});
        theta = new RealParameter(new Double[] {0.5, 0.5});
        nodeMath = new OUNodeMath();
        colorAssignments = new IntegerParameter();
        colorAssignments = new IntegerParameter(new Integer[]{ 0 });
        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta,
                "sigmasq", sigmasq, "covariance", covariance,
                "popVar", popVar, "popCov", popCov,
                "dAlpha", dAlpha,
                "root", rootValues,
                "optNr", 1, "optAssign", colorAssignments);

        pcm = new DOUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(-15.7588578266153, logP, EPSILON);
    }
}
