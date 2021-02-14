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
        treeStr = "(((t3:1.209461463,t4:0.0):0.6659547705,(t9:0.841016425,t10:0.841016425):1.034399809):1.561956365,(((t2:1.602817551,(t5:1.164343725,t6:1.164343725):0.4384738261):0.6643605462,t1:2.267178098):0.7120187616,(t7:1.115655119,t8:1.115655119):1.86354174):0.458175739);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "t3 t4 t9 t10 t2 t5 t6 t1 t7 t8";

        // continuous traits
        nTraits = 2;
        contTraitData = Arrays.asList(
                -0.508635648090195, -0.279620384230921,
                -0.890439381916328, -0.0389259110017274,
                -1.60403656615896, -0.152064627988708,
                -1.78355663933205, 0.0876041569681547,
                -4.55402786822282, 1.72033772810441,
                -1.866891946151, 0.499459863055851,
                -2.78627660312175, 3.17287179990728,
                -7.70315227205359, 3.79277958924578,
                -4.23115951625136, 1.99156992779389,
                -4.99358298983542, 4.48463339235722
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {479.701057745291, 860.744522346309});
        sigmasq = new RealParameter(new Double[] {0.000623020788861537, 7.8185676755921});
        correlation = new RealParameter(new Double[] {0.0440561576525885});
        sigmaMatrix.initByName("sigmasq", sigmasq, "covariance", correlation, "trait", morphData);
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "tree", tree, "rootValues", rootValues);

        alpha = new RealParameter(new Double[] {0.903659679354104, 1.61416481463754, 1.61416481463754, 2.88331082048443});
        theta = new RealParameter(new Double[] {-525940.712324216, 294437.489034769});
        optNr = 1;
        optAssignment = new IntegerParameter(new Integer[] {0});
        modelParameter.initByName("alpha",alpha,"theta", theta, "optNr", optNr, "optAssign", optAssignment, "trait", morphData);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-32.1041358768315, logP, EPSILON);
    }

    @Test
    public void testOULkBDTreeFourTraitsOneOpt() {
        treeStr = "((t2:0.02581382964,((t4:1.364083207,(t6:0.6671516053,(t20:0.01976781459,t21:0.01976781459):0.6473837907):0.6969316013):1.283955909,(t3:0.07558684674,(((((t10:0.1304039739,t11:0.1304039739):0.1802475361,t9:0.3106515099):0.1221429338,((t16:0.07355997008,t17:0.07355997008):0.2046965799,(t14:0.09014646219,t15:0.09014646219):0.1881100878):0.1545378937):0.6837427973,(t5:0.7769685163,((t22:0.017042341,t23:0.017042341):0.4447813857,t7:0.4618237267):0.3151447896):0.3395687247):0.3706573897,(((t18:0.05386849092,t19:0.05386849092):0.2697350422,t8:0.3236035331):0.3127811361,(t12:0.1125485978,t13:0.1125485978):0.5238360714):0.8508099614):0.2945962084):0.8662482765):0.4629658256):0.3514405254,t1:2.468094697);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "t2 t4 t6 t20 t21 t3 t10 t11 t9 t16 t17 t14 t15 t5 t22 t23 t7 t18 t19 t8 t12 t13 t1";

        nTraits = 4;
        contTraitData = Arrays.asList(
                3.00109889967356, -3.70196585575934, 0.944415773151817, -0.627371730519694, 0.705922495817008, -1.1961514290321, -0.365416126164454, 1.22914379964039, 1.06685997301004, 0.949851822464917, -0.895414404665061, -0.26758551373351, 1.47654830127778, 2.51443762671007, 1.72141648936453, -0.288373970192673, 1.79416694284819, 0.999849573611582, 0.913193537813977, 0.223678095393105, 1.63235017771322, -0.845239910823922, -0.0732156266202136, 1.27857592689655, 2.198067603982, 0.677269413717803, 2.49997616153117, 0.263528172162342, 2.00797736357115, -0.583634382464634, 1.34485193822117, 0.334336505660026, 0.837708188055099, -0.439146360303946, 0.314773991228432, -0.261959476735411, 1.70337910382885, 0.584880926753183, 1.85299825183719, 0.860411758912654, 1.03367819169608, 1.25602088839985, 1.18396608579702, 0.920763369782958, 1.57533121713144, 0.298127746103833, 1.30226112912169, 1.56003548885479, 1.34698269753967, 0.767755917031157, 1.00742610052625, 1.65586537515969, 0.736154286398392, 0.848625230273237, 0.653184853339044, 1.29778451177674, 0.377605844641037, 1.77141491316137, -0.152612282629036, 1.86094308483523, 0.713991982767549, 1.4369192590053, -0.154903770005904, 2.11000918063577, 0.524306996575646, 1.52046987619938, 0.0542900235268924, 1.41986855404264, -0.702548110663247, -0.259238914058792, -0.812427627975532, 2.20822905806291, 0.0906773308993026, -1.2496800365324, -0.221056647836304, 2.67862036254945, 0.167572407859486, -1.36419063224393, -0.671092226830047, 2.01909630410049, 1.15780387828872, -0.691169990648929, 2.93916737537487, 2.18803123180906, 0.85641477038843, -1.14669524555737, 1.17637853252482, 3.67378303649472, -0.489514755688319, -1.60280945083744, -0.820563354673034, 3.92823992184346
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {14.9423638086769, -12.8322451408723, -2.12734808227883, 3.0060006335114});
        sigmasq = new RealParameter(new Double[] {1.91889428979701, 6.69847612212845, 4.8862486730982, 1.80310998971464});
        correlation = new RealParameter(new Double[] {-1.69303813349868, 1.31414017941805, 0.768091892963974, 2.42053441234052, -1.25999352981913, -1.26372726667469});
        sigmaMatrix.initByName("sigmasq", sigmasq, "covariance", correlation, "trait", morphData);
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "tree", tree, "rootValues", rootValues);

        alpha = new RealParameter(new Double[] {4.81019687857697, -0.572364908034288, -0.991173642257948, 2.14994820449264, -0.572364908034288, 2.5803246133496, 0.198641987350996, 0.348045384424873, -0.991173642257948, 0.198641987350996, 2.27914049120142, -0.882649865495069, 2.14994820449264, 0.348045384424873, -0.882649865495069, 1.36599586002177});
        theta = new RealParameter(new Double[] {-1.0121752938941, -1.02906717423969, 1.29288300684234, 5.62202767021494});
        optNr = 1;
        optAssignment = new IntegerParameter(new Integer[] {0});
        modelParameter.initByName("alpha",alpha,"theta", theta, "optNr", optNr, "optAssign", optAssignment, "trait", morphData);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-80.5541960395741, logP, EPSILON);
    }

    @Test
    public void testOULkBigTreeThreeTraitsOneOpt(){
        treeStr = "(((((t7:1.170584527,t8:0.06085260709):0.4606487947,t3:0.2614399944):0.579264194,((t4:0.2273832628,((t20:0.2589859644,t21:0.2589859644):1.026733719,t5:1.285719683):0.2750009417):0.09723438815,(((t16:0.4159281643,((t23:0.08051989155,t24:0.08051989155):0.1811578944,t19:0.261677786):0.1542503784):0.3199676779,t12:0.7358958422):0.0005990123775,((t27:0.04888641467,t28:0.04888641467):0.08434037086,t22:0.1332267855):0.603268069):0.9214601586):0.5525425024):0.7579399341,(((((t17:0.3392692305,t18:0.3781764464):0.2709967043,t13:0.6491731507):0.6341627369,(t10:0.8673868947,t11:0.8673868947):0.415948993):0.04839422779,(t6:0.09790727679,(((t14:0.4937862752,t15:0.4937862752):0.03567073945,(t25:0.06689877597,t26:0.06689877597):0.4625582386):0.6021536888,t9:0.8724904128):0.05862586847):0.1414935436):1.256344292,t2:0.2493569227):0.3803630421):0.4341342593,t1:1.156964806);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "t7 t8 t3 t4 t20 t21 t5 t16 t23 t24 t19 t12 t27 t28 t22 t17 t18 t13 t10 t11 t6 t14 t15 t25 t26 t9 t2 t1";

        nTraits = 3;
        contTraitData = Arrays.asList(
                -5.74403666827106, 2.80304937034122, 10.9999817767616, -2.71239841963195, -2.84736645587032, 5.12531545684057, -0.315040453213165, -2.09685157852721, 6.06329301501277, 3.31067677454528, -0.257506250894251, 1.86909015582891, 2.84218731569332, 0.180929066250068, 0.560340536955865, 4.05713487489884, -3.7055443245073, 0.96124618636242, 4.02902233141467, -0.84653926402611, -0.767543506786659, 2.94185215708918, 0.667505918426113, 1.08948513995746, 3.95628963235866, -1.67002916624278, 2.07966137491632, 3.06801288272253, 0.0491340158808997, -0.704703177520037, 2.08915116689399, -2.94666689206912, -1.91539789647769, 2.34474105661159, -0.445818710975407, 2.27917904044858, 2.34134566362511, -0.803642945072774, 0.898920200231013, 2.00677710914333, 2.87810120222232, 0.957445058016412, 2.85851949036752, -0.0844808717397866, 2.76501536913025, 0.968747289084197, 0.701488225639843, 0.697283801738745, 0.0702111605185246, 3.03890576591112, 3.39315988878406, -1.12430942570044, 2.79296064334025, 1.96411215255727, 2.46402739738304, 2.92408958420461, 3.06338712402313, -4.13207227814121, -1.05763196184776, 6.16594708521046, -0.982699854187152, 0.364996941022436, 4.48597895018709, -1.05401942634444, -0.620201908853779, -1.14081385556703, 2.83838409019394, 3.9247861405911, -1.60148745583127, 2.69561412842048, -0.0646631845267889, 2.56487872024003, 0.41350963649962, -1.41993924117363, 1.78466117896427, 1.01556754554521, -0.914895698908689, -1.23127539456333, -3.37333875638925, -0.373858742421192, 4.16052949284078, -1.79225290755084, 0.670041293952786, 7.16605998753728
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {-1.90691945660393, -0.715869367111241, 1.41764335260452});
        sigmasq = new RealParameter(new Double[] {1.6503513882414, 2.24538055884218, 2.08463965942004});
        correlation = new RealParameter(new Double[] {1.0760202029967, -0.680805700219556, 0.837864728261558});
        sigmaesq = new RealParameter(new Double[] {1.17477955284909, 1.22777558870427, 1.39941466666442});
        RealParameter popCov = new RealParameter(new Double[] {0.0402095219320196, -0.150691881961759, 0.130349880732245});
        sigmaMatrix.initByName("sigmasq", sigmasq, "covariance", correlation, "upperMatrix", true, "trait", morphData);
        sigmaEMatrix.initByName("sigmasq", sigmaesq, "covariance", popCov, "upperMatrix", true, "trait", morphData);
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "popMatrix", sigmaEMatrix, "tree", tree, "rootValues", rootValues);

        alpha = new RealParameter(new Double[] {1.12803383714597, 1.40264158623706, 1.144513851248, 0.964737659784515, 1.5083252410321, 0.912498590336744, 0.931791823825421, 1.50189246859842, 0.775617711727039});
        theta = new RealParameter(new Double[] {1.05924392826747, 0.892577097362887, 0.800735679276209});
        optNr = 1;
        optAssignment = new IntegerParameter(new Integer[] {0});
        modelParameter.initByName("alpha",alpha,"theta", theta, "optNr", optNr, "optAssign", optAssignment, "trait", morphData);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-177.637641073355, logP, EPSILON);
    }

}
