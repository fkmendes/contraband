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

public class OULikelihoodTest {
    private final static double EPSILON = 1e-6;
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> contTraitData;
    private RealParameter correlation;
    private RealParameter rootValues;
    private RealParameter sigmasq;
    private RealParameter sigmaesq;
    private RealParameter alpha;
    private RealParameter theta;
    private Integer optNr;
    private IntegerParameter optAssignment;
    private final OUModelParameter modelParameter = new OUModelParameter();
    private final MorphologicalData morphData = new MorphologicalData();
    private final GeneralNodeMath nodeMath = new GeneralNodeMath();
    private final SigmaMatrix sigmaMatrix = new SigmaMatrix();
    private final SigmaMatrix sigmaEMatrix = new SigmaMatrix();
    private final KeyRealParameter contTrait = new KeyRealParameter();
    private final IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
    private final RealParameter colorValues = new RealParameter(new Double[] {1.0});
    private final RateCategoryClockModel lsc = new RateCategoryClockModel();

    @Test
    public void testOULkSmallTreeTwoTraitsOneOpt() {
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

    @Test
    public void testOULkSmallTreeThreeTraitsOneOpt() {
        treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "sp1 sp2 sp3 sp4";

        // continuous traits
        nTraits = 3;
        contTraitData = Arrays.asList(
                2.90115170364898, 2.14268145872343, 3.2824264362529,
                2.62770319717086, 1.97031352144559, 4.28940778474205,
                3.96814434248716, 3.19423412470143, 3.49323574837925,
                8.88546166300816, 9.77485338099449, -3.65830173123419
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        sigmasq = new RealParameter(new Double[] {0.35, 0.2, 0.8});
        correlation = new RealParameter(new Double[] {0.06, 0.06, 0.4});
        sigmaMatrix.initByName("sigmasq", sigmasq, "covariance", correlation, "upperMatrix", true, "trait", morphData);
        rootValues = new RealParameter(new Double[] {2.2, 1.3, 0.5});
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "tree", tree, "rootValues", rootValues);

        alpha = new RealParameter(new Double[] {1.8, -0.9, 1.2, -0.9, 1.0, 1.4, 1.2, 1.4, 1.6});
        theta = new RealParameter(new Double[] {5.0, 4.8, 1.2});
        optNr = 1;
        optAssignment = new IntegerParameter(new Integer[] {0});
        modelParameter.initByName("alpha",alpha,"theta", theta, "optNr", optNr, "optAssign", optAssignment, "trait", morphData);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-260.827821834006, logP, EPSILON);
    }

    @Test
    public void testOULkSmallTreeTwoTraitsTwoOpt(){
        treeStr = "(((sp1:2,sp2:2):3,sp3:5):1,sp4:6);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "sp1 sp2 sp3 sp4";

        // continuous traits
        nTraits = 2;
        contTraitData  = Arrays.asList(
                -0.852088501326528, -0.724973069661161,
                0.0490989886775859, 0.478463675183718,
                0.477307495634625, -0.0761208848808844,
                0.826579895044279, 0.122281759527914
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // OU model parameters
        sigmasq = new RealParameter(new Double[] {1.0, 1.0});
        correlation = new RealParameter(new Double[] {0.6});
        sigmaMatrix.initByName("sigmasq", sigmasq, "covariance", correlation, "upperMatrix", true, "trait", morphData);
        rootValues = new RealParameter(new Double[] {0.2, 1.3});
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "tree", tree, "rootValues", rootValues);

        alpha = new RealParameter(new Double[] {1.0, 1.2, 1.2, 1.0});
        theta = new RealParameter(new Double[] {4.2, 3.8, 1.0, 2.0});
        optNr = 2;
        optAssignment = new IntegerParameter(new Integer[]{ 1, 1, 0, 0, 0, 0, 0});
        modelParameter.initByName("alpha",alpha,"theta", theta, "optNr", optNr, "optAssign", optAssignment, "trait", morphData);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-95.4784626512006, logP, EPSILON);
    }

    @Test
    public void testOULkSATreeTwoTraitsOneOpt(){

    }

}
