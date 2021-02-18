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
    public void testJOULkSmallTreeTwoTraitsOneOptTwoJumps() {
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
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "popMatrix", sigmaEMatrix, "tree", tree, "rootValues", rootValues);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-7.69394705414663, logP, EPSILON);
    }

    @Test
    public void testJOULkLargeTreeThreeTraitsOneOptThreeJumps (){
        treeStr = "((((t12:0.760334619,(t17:0.1354838491,(t20:0.04258159395,t21:0.04258159395):0.09290225512):0.6248507699):0.1501356711,(t9:0.8644970359,(t15:0.1632534095,t16:0.1632534095):0.7012436264):0.04597325419):1.082500367,(t4:0.07133049472,(t7:0.9277357945,t8:0.2883032296):0.4259908727):0.6392439899):1.09942315,(t1:2.959344055,((t13:0.7086026648,t14:0.7086026648):1.072697974,(((t2:0.595474196,t3:0.1109147523):0.1540832003,(t5:0.6520365485,(t6:1.02081611,(t10:0.8045108475,t11:0.8045108475):0.4490546309):0.09812952038):0.1880926274):0.1385491448,(t18:0.05460515063,t19:0.05460515063):1.62373162):0.1029638676):1.178043416):0.1330497523);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "t12 t17 t20 t21 t9 t15 t16 t4 t7 t8 t1 t13 t14 t2 t3 t5 t6 t10 t11 t18 t19";

        nTraits = 3;
        contTraitData = Arrays.asList(
                -0.943581281640293, -0.481731693118047, -0.770020592157489, -1.96815432492033, -0.316578637149189, 1.96293182978666, -0.941014060242354, 0.0197842754383456, 2.82757959444878, -0.815456182012012, -0.0641217507082689, 2.78363317459176, -0.444501575156234, 0.258016621807697, -0.304731942952906, 0.186424440534394, 0.233602570379118, -0.112867020723581, -0.990840681601696, 0.156466577818302, 0.903473276231531, 0.111261673272528, 0.549741570702429, 0.452403717272718, -0.287680334916718, 0.566869411079293, 1.04502265360728, -0.464559560532412, 0.38241148215889, 0.774590130690366, -0.655471572889876, 0.450582026890006, 1.2271752135018, -0.14896092655346, 0.0980894351987711, -2.30308374280533, -0.436942473744585, 0.593134147695387, 2.00629833548485, 0.326844662990295, 0.493149271034999, 0.0288925535767135, 0.0870931645759436, 0.724421143245464, 0.844840081258482, 0.755837904922845, 0.849462895085013, -0.190775160210953, 0.0381135613925536, 0.696489622440541, -0.689504730941722, 0.652399311714125, 0.583230324495994, -1.28586558702587, 0.444798330841347, 0.50054076113854, -1.10039140013272, 0.1430044914221, 0.63585861819614, -0.419374943995031, 0.389587679086552, 0.577264264586204, -1.21436211611273
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {0.3, 0.4, 0.2});
        sigmasq = new RealParameter(new Double[] {1.0, 0.2, 2.0});
        covariance = new RealParameter(new Double[] {0.6, 0.3, 0.5});
        sigmaMatrix.initByName("sigmasq", sigmasq, "covariance", covariance, "upperMatrix", true, "trait", morphData);

        alpha = new RealParameter(new Double[] {2.0, -1.0, 1.0, -1.0, 3.0, 0.5, 1.0, 0.5, 1.5});
        theta = new RealParameter(new Double[] {0.4, 0.6, -0.5});
        jumpVar = new RealParameter(new Double[] {0.4, 0.1, 0.2});
        jumpCov = new RealParameter(new Double[] {0.6, 0.3, 0.5});
        sigmaJMatrix.initByName("sigmasq", jumpVar, "covariance", jumpCov,  "upperMatrix", true,"trait", morphData);
        jumpMean = new RealParameter(new Double[] {0.1, 0.2, -0.3});
        jumpIndicators = new IntegerParameter(new Integer[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
        optNr = 1;
        optAssignment = new IntegerParameter(new Integer[] {0});
        modelParameter.initByName("alpha",alpha,"theta", theta, "optNr", optNr, "optAssign", optAssignment, "trait", morphData, "jumpVCVMat", sigmaJMatrix, "jumpMean", jumpMean, "jumpIndicator", jumpIndicators);
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "tree", tree, "rootValues", rootValues);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-38.7133704437508, logP, EPSILON);
    }

}
