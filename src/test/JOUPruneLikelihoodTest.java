package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.prunelikelihood.OUNodeMath;
import contraband.prunelikelihood.JOUPruneLikelihood;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class JOUPruneLikelihoodTest {
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
    RealParameter jumpVar;
    RealParameter jumpCov;
    RealParameter alpha;
    RealParameter theta;
    RealParameter jumpMean;
    IntegerParameter jumpIndicators;
    IntegerParameter colorAssignments;
    JOUPruneLikelihood pcm;

    /*
     * (1) Ultrametric tree with 3 taxa, 2 continuous traits, 1 optima
     */
    @Test
    public void testOUPruneLkSmallTreeTwoTraitsOneOpt() {
        treeStr = "((A:3.0,B:3.0):4.0,C:7.0);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "A B C";

        nTraits = 2;
        oneTraitValues = Arrays.asList(
                0.143769382168306, -0.19301903320227,
                -0.661874557048517, -0.360841881116811,
                -0.216063301731694, 1.28966779981295
        );
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {0.3, 0.4});
        sigmasq = new RealParameter(new Double[] {1.0, 1.0});
        covariance = new RealParameter(new Double[] {0.6});
        popVar = new RealParameter(new Double[] {1.0, 1.0});
        popCov = new RealParameter(new Double[] {0.4});
        alpha = new RealParameter(new Double[] {2.0, 1.0, 1.0, 2.0});
        theta = new RealParameter(new Double[] {0.5, 0.5});
        jumpVar = new RealParameter(new Double[] {1.0, 1.0});
        jumpCov = new RealParameter(new Double[] {0.2});
        jumpMean = new RealParameter(new Double[] {0.1, 0.2});
        jumpIndicators = new IntegerParameter(new Integer[]{1, 0, 0, 1, 0});
        nodeMath = new OUNodeMath();
        colorAssignments = new IntegerParameter();
        colorAssignments = new IntegerParameter(new Integer[]{ 0 });
        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta,
                "sigmasq", sigmasq, "covariance", covariance,
                "popVar", popVar, "popCov", popCov,
                "root", rootValues,
                "jumpVar", jumpVar, "jumpCov", jumpCov, "jumpMean", jumpMean, "jump", jumpIndicators,
                "optNr", 1, "optAssign", colorAssignments);

        pcm = new JOUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(-7.69394705414663, logP, EPSILON);
    }
}
