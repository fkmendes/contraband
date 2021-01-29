package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.NodeMath;
import contraband.prunelikelihood.LiabilityLikelihood;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class LiabilityLikelihoodTest {

    private final static double EPSILON = 1e-6;
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> data;
    private RealParameter inverseMatrix;
    private RealParameter sigmasq;
    private final NodeMath nodeMath = new NodeMath();
    private final KeyRealParameter traitValues = new KeyRealParameter();
    private final IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
    private final RealParameter colorValues = new RealParameter(new Double[] {1.0});
    private final RateCategoryClockModel lsc = new RateCategoryClockModel();

    @Test
    public void testLiabilityLikelihoodInverseMatrix() {
        // tree with sampled ancestors
        treeStr = "(((t3:1.209461463,t4:0.0):0.6659547705,(t9:0.841016425,t10:0.841016425):1.034399809):1.561956365,(((t2:1.602817551,(t5:1.164343725,t6:1.164343725):0.4384738261):0.6643605462,t1:2.267178098):0.7120187616,(t7:1.115655119,t8:1.115655119):1.86354174):0.458175739);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "t3 t4 t9 t10 t2 t5 t6 t1 t7 t8";

        // trait values
        nTraits = 2;
        data = Arrays.asList(
                -4.12523015629711, 5.46820550439604,
                -1.62252112839929, 1.33650033921425,
                -2.44641157681447, 2.84921644558551,
                -1.54740433331906, -0.0963333650163727,
                0.881613016929177, -1.79728580710375,
                -0.217262173674136, 1.72609075255175,
                1.62614740980059, -1.7120872848909,
                -2.89454199395993, 4.49117923198017,
                2.45457588141588, -6.28243322556914,
                2.38497975048708, -5.5445431329235
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // trait evolutionary rate matrix -> inverse matrix
        sigmasq = new RealParameter(new Double[]{1.0});
        inverseMatrix = new RealParameter(new Double[] {7.96207107912967, 4.20539442730983, 4.20539442730983, 2.42683976729003});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "inverseMatrix", inverseMatrix);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // likelihood
        LiabilityLikelihood PCM1 = new LiabilityLikelihood();
        PCM1.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lnLk = PCM1.calculateLogP();
        Assert.assertEquals(-33.6679664296583, lnLk, EPSILON);
    }
}
