package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeTrait;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.GeneralNodeMath;
import contraband.prunelikelihood.MorphologicalData;
import contraband.prunelikelihood.OULikelihood;
import contraband.prunelikelihood.OUModelParameter;
import contraband.prunelikelihood.SigmaMatrix;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

public class OULikelihoodTest {
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
    private RealParameter alpha;
    private RealParameter theta;
    private RealParameter delta;
    private RealParameter lambda;
    private RealParameter inverseMatrix;
    private Integer optNr;
    private IntegerParameter optAssignment;
    private final OUModelParameter modelParameter = new OUModelParameter();
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
        treeStr = "((A:3.0058179,B:3.0058179):4.350951,C:7.3567689);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "A B C";

        // continuous traits
        nTraits = 2;
        contTraitData = Arrays.asList(
                0.131394584822684, -0.19269144948362,
                -0.65643500381027,-0.373259036803319,
                -0.17399894787897, 1.27056761078824
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // trait rate matrix and population variance
        sigmasq = new RealParameter(new Double[] {1.36, 1.0});
        correlation = new RealParameter(new Double[] {0.6});
        sigmaesq = new RealParameter(new Double[] {1.09, 1.0});
        RealParameter popCov = new RealParameter(new Double[] {0.3});
        sigmaMatrix.initByName("sigmasq", sigmasq, "covariance", correlation, "trait", morphData);
        sigmaEMatrix.initByName("sigmasq", sigmaesq, "covariance", popCov, "trait", morphData);
        rootValues = new RealParameter(new Double[] {0.2, 1.3});
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "popMatrix", sigmaEMatrix, "tree", tree, "rootValues", rootValues);

        alpha = new RealParameter(new Double[] {2.0, 0.0, 0.0, 2.0});
        theta = new RealParameter(new Double[] {0.5, 0.5});
        optNr = 1;
        optAssignment = new IntegerParameter(new Integer[] {0});
        modelParameter.initByName("alpha",alpha,"theta", theta, "optNr", optNr, "optAssign", optAssignment, "trait", morphData);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-7.58111239313721, logP, EPSILON);
    }


}
