package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import contraband.prunelikelihood.OUPruneLikelihood;
import org.apache.commons.math3.linear.*;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.prunelikelihood.OUNodeMath;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

public class OUNodeMathTest {
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
    OUPruneLikelihood pcm;
    IntegerParameter colorAssignments;

    @Test
    public void testAlphaDecomposition() {
        treeStr = "((A:3.0058179,B:3.0058179):4.350951,C:7.3567689);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "A B C";

        nTraits = 2;
        oneTraitValues = Arrays.asList(
                4.53371989048144, 4.39037816144705,
                4.64369804601356, 4.50915779102972,
                6.22349219457299, 6.60199054481913
        );
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {0.2, 1.3});
        sigma = new RealParameter(new Double[] {1.36, 0.6, 1.0});
        sigmae = new RealParameter(new Double[] {1.0, 0.3, 1.0});
        alpha = new RealParameter(new Double[] {1.0, 3.0, 2.0, 4.0});
        theta = new RealParameter(new Double[] {0.5, 0.5});

        colorAssignments = new IntegerParameter();
        colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0 });

        nodeMath = new OUNodeMath();
        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta,
                "sigma", sigma, "sigmae", sigmae,
                "root", rootValues,"optNr", 1, "optAssign", colorAssignments,
                "upperMatrix", false);

        pcm = new OUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);


        double logP = pcm.calculateLogP();
        Assert.assertEquals(-56.447173900514, logP, EPSILON);


/*
        nodeMath.getLMatForNode(0);
        nodeMath.getLMatForNode(1);
        nodeMath.getLMatForNode(2);
        nodeMath.getLMatForNode(3);
        nodeMath.getLMatForNode(4);

        nodeMath.getMVecForNode(0);
        nodeMath.getMVecForNode(1);
        nodeMath.getMVecForNode(2);
        nodeMath.getMVecForNode(3);
        nodeMath.getMVecForNode(4);

        nodeMath.getRForNode(0);
        nodeMath.getRForNode(1);
        nodeMath.getRForNode(2);
        nodeMath.getRForNode(3);
        nodeMath.getRForNode(4);

        nodeMath.getAMatForNode(0);
        nodeMath.getAMatForNode(1);
        nodeMath.getAMatForNode(2);
        nodeMath.getAMatForNode(3);
        nodeMath.getAMatForNode(4);

        nodeMath.getCMatForNode(0);
        nodeMath.getCMatForNode(1);
        nodeMath.getCMatForNode(2);
        nodeMath.getCMatForNode(3);
        nodeMath.getCMatForNode(4);
*/
    }

}
