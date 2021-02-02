package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.StandardData;
import beast.evolution.datatype.UserDataType;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.NodeMath;
import contraband.prunelikelihood.BinaryDiscreteTraits;
import contraband.prunelikelihood.LiabilityLikelihood;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;

import java.util.ArrayList;
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
    private RealParameter correlation;
    private RealParameter rootValues;
    private final NodeMath nodeMath = new NodeMath();
    private final KeyRealParameter traitValues = new KeyRealParameter();
    private final IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
    private final RealParameter colorValues = new RealParameter(new Double[] {1.0});
    private final RateCategoryClockModel lsc = new RateCategoryClockModel();
/*
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
*/
    @Test
    public void testContBinaryTraitsLikelihood() {
        // tree with sampled ancestors
        treeStr = "(((t3:1.275476423,(t4:0.5941849179,(t5:0.3540775541,t6:0.3540775541):0.2401073638):0.6812915054):0.02631446592,t2:1.301790889):1.054321576,t1:2.356112466);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "t3 t4 t5 t6 t2 t1";

        // trait values
        nTraits = 1;
        data = Arrays.asList(
                0.82814124000102, -1.09500909925769, 1.15728364543577, 1.52672957333196, 1.69564041037902, -2.7406878251688
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // each species has four binary discrete traits
        List<Sequence> sequenceList = new ArrayList<>(6);
        Sequence sp1 = new Sequence("t1", "1");
        Sequence sp2 = new Sequence("t2", "0");
        Sequence sp3 = new Sequence("t3", "0");
        Sequence sp4 = new Sequence("t4", "1");
        Sequence sp5 = new Sequence("t5", "0");
        Sequence sp6 = new Sequence("t6", "0");
        sequenceList.add(0, sp1);
        sequenceList.add(1, sp2);
        sequenceList.add(2, sp3);
        sequenceList.add(3, sp4);
        sequenceList.add(4, sp5);
        sequenceList.add(5, sp6);

        UserDataType userDataType1 = new UserDataType();
        userDataType1.initByName("characterName", "ch1","codeMap", "0=0, 1=1", "states", 2, "value", "0 Present, 1 Absent");

        List<UserDataType> charStateLabels= new ArrayList<>(1);
        charStateLabels.add(0, userDataType1);

        StandardData standardData = new StandardData();
        standardData.initByName("charstatelabels", charStateLabels);

        Alignment data = new Alignment();
        data.initByName("sequence", sequenceList, "userDataType", standardData);

        RealParameter liabilities  = new RealParameter(new Double[] {0.611885645256119, -1.21656277758523, -1.27268374351877, 0.426319816642301, -0.439892339448248, -0.302559487770471});
        BinaryDiscreteTraits discreteTraits = new BinaryDiscreteTraits();
        discreteTraits.initByName("liability", liabilities, "data", data, "tree", tree);

        discreteTraits.calculateLogP();
        double lnLkDist = discreteTraits.getLogP();
        Assert.assertEquals(0.0, lnLkDist, EPSILON);

        // trait evolutionary rate matrix -> inverse matrix
        sigmasq = new RealParameter(new Double[]{1.82570004166763, 1.0});
        correlation = new RealParameter(new Double[] {-0.800315997780182});
        rootValues = new RealParameter(new Double[]{0.284038569053525, -0.433332287053399});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "binaryDiscreteTraits", discreteTraits, "correlation", correlation, "rootValues", rootValues);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // likelihood
        LiabilityLikelihood PCM1 = new LiabilityLikelihood();
        PCM1.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "binaryDiscreteTraits", discreteTraits, "branchRateModel", lsc);
        double lnLk = PCM1.calculateLogP();
        Assert.assertEquals(-16.7994076650072, lnLk, EPSILON);
    }
}
