package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.GeneralNodeMath;
import contraband.prunelikelihood.*;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class JOULikelihoodTest {
    private final static double EPSILON = 1e-6;
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> contTraitData;
    private RealParameter covariance;
    private RealParameter rootValues;
    private RealParameter sigmasq;
    private RealParameter popVar;
    private RealParameter popCov;
    private RealParameter jumpVar;
    private RealParameter jumpCov;
    private RealParameter jumpMean;
    private IntegerParameter jumpIndicators;
    private RealParameter alpha;
    private RealParameter theta;
    private Integer optNr;
    private IntegerParameter optAssignment;
    private final OUModelParameter modelParameter = new OUModelParameter();
    private final MorphologicalData morphData = new MorphologicalData();
    private final GeneralNodeMath nodeMath = new GeneralNodeMath();
    private final SigmaMatrix sigmaMatrix = new SigmaMatrix();
    private final SigmaMatrix sigmaEMatrix = new SigmaMatrix();
    private final SigmaMatrix sigmaJMatrix = new SigmaMatrix();
    private final KeyRealParameter contTrait = new KeyRealParameter();
    private final IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
    private final RealParameter colorValues = new RealParameter(new Double[] {1.0});
    private final RateCategoryClockModel lsc = new RateCategoryClockModel();

    @Test
    public void testJOUPruneLkSmallTreeTwoTraitsOneOptTwoJumps() {
        treeStr = "((A:3.0,B:3.0):4.0,C:7.0);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "A B C";

        nTraits = 2;
        contTraitData = Arrays.asList(
                0.143769382168306, -0.19301903320227,
                -0.661874557048517, -0.360841881116811,
                -0.216063301731694, 1.28966779981295
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {0.3, 0.4});
        sigmasq = new RealParameter(new Double[] {1.0, 1.0});
        covariance = new RealParameter(new Double[] {0.6});
        sigmaMatrix.initByName("sigmasq", sigmasq, "covariance", covariance, "upperMatrix", true, "trait", morphData);

        popVar = new RealParameter(new Double[] {1.0, 1.0});
        popCov = new RealParameter(new Double[] {0.4});
        sigmaEMatrix.initByName("sigmasq", popVar, "covariance", popCov,  "upperMatrix", true,"trait", morphData);

        jumpVar = new RealParameter(new Double[] {1.0, 1.0});
        jumpCov = new RealParameter(new Double[] {0.2});
        sigmaJMatrix.initByName("sigmasq", jumpVar, "covariance", jumpCov,  "upperMatrix", true,"trait", morphData);

        jumpMean = new RealParameter(new Double[] {0.1, 0.2});
        jumpIndicators = new IntegerParameter(new Integer[]{1, 0, 0, 1, 0});

        alpha = new RealParameter(new Double[] {2.0, 1.0, 1.0, 2.0});
        theta = new RealParameter(new Double[] {0.5, 0.5});
        optNr = 1;
        optAssignment = new IntegerParameter(new Integer[] {0});
        modelParameter.initByName("alpha",alpha,"theta", theta, "optNr", optNr, "optAssign", optAssignment, "trait", morphData, "jumpVCVMat", sigmaJMatrix, "jumpMean", jumpMean, "jumpIndicator", jumpIndicators);
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "popMatrix", sigmaEMatrix, "jumpMatrix", sigmaJMatrix, "tree", tree, "rootValues", rootValues);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-7.69394705414663, logP, EPSILON);
    }
}
