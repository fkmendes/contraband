package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import org.junit.Assert;
import org.junit.Test;
import contraband.math.NodeMath;
import contraband.prunelikelihood.BMPruneLikelihood;
import java.util.Arrays;
import java.util.List;

/*
 * This class contains unit tests for BMPruneLikelihood.class
 * (without using shrinkage method)
 */

public class BMPruneLikelihoodTest {

    private final static double EPSILON = 1e-6;
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> data;
    private RealParameter correlation;
    private RealParameter sigmasq;
    private final NodeMath nodeMath = new NodeMath();
    private final RealParameter traitValues = new RealParameter();
    private final IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
    private final RealParameter colorValues = new RealParameter(new Double[] {1.0});
    private final RateCategoryClockModel lsc = new RateCategoryClockModel();

    /*
     * tree with 3 species and each species has 1 continuous trait
     */
    @Test
    public void testBMPruneLikelihood3Species1Trait() {
        // tree
        treeStr = "((sp1:1.0, sp2:1.0):1.0,sp3:2.0);";
        spNames = "sp1 sp2 sp3";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 1;
        data = Arrays.asList(
                0.3022866,
                0.4237119,
                2.5298873
        );
        traitValues.initByName("value", data, "keys", spNames);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // BM model parameters
        sigmasq = new RealParameter(new Double[]{0.449638667542283});
        correlation = new RealParameter(new Double[]{0.0});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation);

        // prune likelihood
        BMPruneLikelihood fpcm1 = new BMPruneLikelihood();
        fpcm1.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lik1 = fpcm1.calculateLogP();

        Assert.assertEquals(-3.95372886237383, lik1, EPSILON);
    }

    /*
     * tree with 10 species and each species has 2 continuous traits
     */
    @Test
    public void testBMPruneLikelihood10Species2Traits() {
        // tree
        treeStr = "(((t3:1.209461463,t4:1.209461463):0.6659547705,(t9:0.841016425,t10:0.841016425):1.034399809):1.561956365,(((t2:1.602817551,(t5:1.164343725,t6:1.164343725):0.4384738261):0.6643605462,t1:2.267178098):0.7120187616,(t7:1.115655119,t8:1.115655119):1.86354174):0.458175739);";
        spNames = "t3 t4 t9 t10 t2 t5 t6 t1 t7 t8";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 2;
        data =  Arrays.asList(
                0.326278727608277, -3.22668212260941,
                1.8164550628074, -1.71183724870188,
                -0.370085503473201, 1.81925405275285,
                0.665116986641999, -0.428821390843389,
                1.17377224776421, 4.22298205455098,
                3.59440970719762, 1.51483058860744,
                3.38137444987329, 3.63674837097173,
                -0.187743059073837, 3.68456953445085,
                -1.64759474375234, -0.743303344769609,
                -2.19534387982435, 1.10602125889508
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // BM model parameters
        sigmasq = new RealParameter(new Double[]{0.93818618860063, 1.74928729364035});
        correlation = new RealParameter(new Double[]{-0.117852694436286});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation);

        // prune likelihood
        BMPruneLikelihood fpcm2 = new BMPruneLikelihood();
        fpcm2.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lik2 = fpcm2.calculateLogP();
        Assert.assertEquals(-39.588461762377, lik2, EPSILON);
    }

    /*
     * tree with 12 species and each species has 4 continuous traits
     */
    @Test
    public void testBMPruneLikelihood12Species4Traits() {
        // tree
        treeStr = "(t1:0.4387394809,((((t9:0.04695960989,t10:0.04695960989):0.3330763616,t3:0.3800359715):0.02593808098,(t4:0.3763515397,((t5:0.3336782345,t6:0.3336782345):0.001793424038,((t11:0.01775546957,t12:0.01775546957):0.2525263158,(t7:0.1384768423,t8:0.1384768423):0.1318049431):0.06518987311):0.04087988122):0.02962251275):0.005880353322,t2:0.2469241568):0.3393625325);";
        spNames = "t1 t9 t10 t3 t4 t5 t6 t11 t12 t7 t8 t2";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 4;
        data =  Arrays.asList(
                0.353663344600893, 0.229894452237773, -1.04988557964476, -1.00377460950264,
                -0.260659762787646, -0.889968417857947, -1.55616727235143, -1.58433907435323,
                -0.143297771543009, -1.01928746070036, -1.09488792059502, -1.15487441742173,
                -0.320198017359097, -0.221714820148973, -0.316592979435508, 0.193491613863128,
                1.58652065854415, 0.866034091552559, 0.607662577851945, 0.858965417527042,
                0.473150703285599, -0.547837172102502, 1.05672142780744, 0.435337697033783,
                1.13940033793728, -0.0992925133953159, 2.11008014995751, 1.67013615353965,
                0.181057545027617, -0.10664389071908, 0.622428141848601, 0.461563348647148,
                0.161848086456716, -0.351658729601069, 0.274713820329526, 0.0131752572461055,
                0.876042348041411, -0.391407502463039, 0.934277091068923, 0.628852350924242,
                1.34896801425551, -0.742998450279645, 1.19164140228254, 0.649028779617329,
                0.735955677879521, -0.974386073214283, 0.120739155255668, 0.00815911328976909
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // BM model parameters
        sigmasq = new RealParameter(new Double[]{0.8048720, 0.8008982, 2.0987084, 1.8666865});
        correlation = new RealParameter(new Double[]{0.249601983702492, 0.621690087985393, 0.557763578318542, 0.253263097787055, 0.465328439817149, 0.942027680951926});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation);

        // prune likelihood
        BMPruneLikelihood fpcm3 = new BMPruneLikelihood();
        fpcm3.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lik3 = fpcm3.calculateLogP();
        Assert.assertEquals(-24.3878360967573, lik3, EPSILON);
    }

    /*
     * tree with 3 species and each species has 1 continuous trait
     * the tree has 1 sampled ancestors, i.e. t2
     */

    @Test
    public void testBMPruneLikelihood3Species1TraitSATree() {
        // tree
        treeStr = "((t1:1.0,t2:0.0):2.0,t3:3.0);";
        spNames = "t1 t2 t3";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 1;
        data = Arrays.asList(
                2.0,3.0,1.0
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // BM model parameters
        correlation =  new RealParameter(new Double[] {
                0.0
        });
        sigmasq = new RealParameter(new Double[] {
                0.6
        });
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // prune likelihood
        BMPruneLikelihood PCM2 = new BMPruneLikelihood();
        PCM2.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lnLk = PCM2.calculateLogP();

        Assert.assertEquals(-4.38645689857906, lnLk, EPSILON);
    }

    /*
     * tree with 10 species and each species has 2 continuous traits
     * the tree has 1 sampled ancestors, i.e. t4
     */

    @Test
    public void testBMPruneLikelihood10Species2TraitsSATree() {
        // tree
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

        // BM model parameters
        correlation = new RealParameter(new Double[] {
                -0.956694237975469
        });
       sigmasq = new RealParameter(new Double[] {
                1.48219484531245,
                4.86284297404418
        });
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        BMPruneLikelihood PCM1 = new BMPruneLikelihood();
        PCM1.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lnLk = PCM1.calculateLogP();

        Assert.assertEquals(-33.6679664296583, lnLk, EPSILON);
    }

    /*
     * tree with 12 species and each species has 3 continuous traits
     * the tree has 1 sampled ancestors, i.e. t4
     */

    @Test
    public void testBMPruneLikelihood12Species3TraitsSATree() {
        // tree
        treeStr = "(t1:0.0162595512,((t2:0.2981734165,t3:0.07231244369):0.320819053,((((t11:0.003069371137,t12:0.003069371137):0.2591040945,t7:0.2621734657):0.1989762127,(t8:0.1021021305,(t9:0.02930302334,t10:0.02930302334):0.07279910718):0.3590475478):0.02879665668,(t4:0.0,(t5:0.2669991903,t6:0.2669991903):0.049200328):0.1737468167):0.3260942874):0.1839593776):0.0;";
        spNames = "t1 t2 t3 t11 t12 t7 t8 t9 t10 t4 t5 t6";
        tree = new TreeParser(treeStr, false, false, true, 0);

        nTraits = 3;
        data = Arrays.asList(
                0.215000043472759, 0.141278178809372, 0.226642733376528,
                0.776521767392364, 0.668704872177164, 0.176856072807917,
                0.274864229944774, 0.806715774880128, -0.762567216544748,
                -0.466686896787782, -0.265669698753016, 0.0870906123165553,
                -0.588015333642268, -0.332967508179201, -0.0827301327675375,
                -0.262434175362964, 0.031825032167221, 0.112420939261318,
                -0.87131974904876, -1.16155653990928, 0.00307828088038135,
                -1.16217914558771, -1.41785772734122, -0.668474831682343,
                -1.43173438937441, -1.99988027886332, -0.768766703759037,
                -0.232384120150305, -1.13606986874854, 0.391361803506349,
                -0.217673646874237, -0.68027906300786, 0.317349178510409,
                -0.271431937473855, -1.65892990236907, 0.0162183133279696
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // BM model parameters
        correlation = new RealParameter(new Double[] {
                0.701335423789108, 0.752711765275858, 0.281971103959096
        });
        sigmasq = new RealParameter(new Double[] {
                0.638576978038351, 1.39965531927086, 1.05567639681806
        });
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // prune likelihood
        BMPruneLikelihood PCM3 = new BMPruneLikelihood();
        PCM3.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lnLk = PCM3.calculateLogP();

        Assert.assertEquals(-9.5414237837198, lnLk, EPSILON);
    }

    /*
     * tree with 11 species and each species has 3 continuous traits
     * the tree has 2 sampled ancestors, i.e. t1 and t5
     */

    @Test
    public void testBMPruneLikelihood11Species3TraitsSATree(){
        // tree
        treeStr = "(((((t3:0.1588129446,t4:0.1588129446):0.01705137773,((t6:0.09203097793,t7:0.09203097793):0.02043978126,t5:0.0):0.06339356319):0.06941873281,(t8:0.03465123152,t9:0.03465123152):0.2106318237):0.4596163207,t1:0.0):0.4342545645,((t10:0.01280757,t11:0.01280757):0.6723641561,t2:0.6851717261):0.4539822143);";
        spNames = "t3 t4 t6 t7 t5 t8 t9 t1 t10 t11 t2" ;
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 3;
        data = Arrays.asList(
                -1.21564787174142, -4.03137227081675, -0.658138349771025,
                -1.1851671823634, -3.60553672704412, -0.0445842979883116,
                -0.694393171280324, -3.03652095560971, 0.521821220721283,
                -0.590751944909478, -3.44783973931706, -0.596693734056411,
                -0.884946580491066, -3.48302742703454, -0.173523045052619,
                -0.639187501019874, -2.23470805810516, -0.0197703417294821,
                -0.705615140755287, -2.06463132261412, 0.136084231018102,
                -0.543347975313145, -1.45735336840457, -0.883187635997939,
                -0.41310125318572, 2.02935389582742, -0.398893370910344,
                -0.578200681462636, 2.04503536970991, -0.596179783025397,
                0.10299350447178, 0.793454794733989, 0.407591907174602
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // BM model parameters
        correlation = new RealParameter(new Double[] {
                0.431864494796763, 0.269799734818774, 0.265743456802657
        });
        sigmasq = new RealParameter(new Double[] {
                0.330573809129692, 2.10107882032566, 1.16247744000489
        });
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // prune likelihood
        BMPruneLikelihood PCM5 = new BMPruneLikelihood();
        PCM5.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lnLk = PCM5.calculateLogP();

        Assert.assertEquals(-19.6436871314934, lnLk, EPSILON);
    }

    /*
     * tree with 12 species and each species has 4 continuous traits
     * the trait rate matrix is populated by an upper matrix, i.e.
     * Sigma <- t(sigma) * sigma
     */

    @Test
    public void testBMPruneLikelihoodWithUpperMatrix12Species4Traits() {
        // tree
        treeStr = "(t1:0.4387394809,((((t9:0.04695960989,t10:0.04695960989):0.3330763616,t3:0.3800359715):0.02593808098,(t4:0.3763515397,((t5:0.3336782345,t6:0.3336782345):0.001793424038,((t11:0.01775546957,t12:0.01775546957):0.2525263158,(t7:0.1384768423,t8:0.1384768423):0.1318049431):0.06518987311):0.04087988122):0.02962251275):0.005880353322,t2:0.2469241568):0.3393625325);";
        spNames = "t1 t9 t10 t3 t4 t5 t6 t11 t12 t7 t8 t2";
        tree = new TreeParser(treeStr, false, false, true, 0);

        //trait values
        nTraits = 4;
        data = Arrays.asList(
                -0.0733746173961214, -0.381897460527792, 3.01043044622191, -0.377856855072646,
                5.68659715338688, -5.45349811184419, -3.77414740185529, 0.57882637837598,
                4.63420888471096, -5.88280081877496, -3.4666212861778, -0.148236089777769,
                -0.318730233572463, -2.70069110253073, -2.53053047827026, 0.862306413264445,
                1.62978685097024, -6.84755589169306, 0.481033646574288, 0.53010180261091,
                0.376724598790307, -0.440758737091378, 2.01600437393992, 2.00282132665632,
                3.00252253755152, -3.3301298265687, -0.367727022453056, 0.657269080077358,
                0.398482000569861, -3.82369481378681, 0.0209910108874771, -0.765003421521001,
                2.22178666226494, -4.61332643056133, -1.12402308648343, -0.857024153578214,
                0.180126949417046, -0.193404866088525, 0.791628971761537, 0.906574999646124,
                0.591666553277025, -1.45443665996485, 1.99527559017226, 0.605770740459275,
                4.94126662089152, -10.453621232694, -1.72287824100249, -0.133159221643173
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // BM model parameters
        sigmasq = new RealParameter(new Double[] {
                2.82578533469758, 4.54723563720239, 3.12137890774279, 1.85982357228878
        });
        correlation = new RealParameter(new Double[] {
                -0.50782453129068, -0.915880932938308, -0.344158561434597,
                0.909007298294455, 0.779078632127494,
                0.385606812313199
        });
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation, "upperMatrix", true);

        // prune likelihood
        BMPruneLikelihood fpcm1 = new BMPruneLikelihood();
        fpcm1.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lnLk = fpcm1.calculateLogP();

        Assert.assertEquals(-86.0358011352751, lnLk, EPSILON);
    }

    /*
     * tree with 13 species and each species has 10 continuous traits
     * the trait rate matrix is populated by an upper matrix, i.e.
     * Sigma <- t(sigma) * sigma
     */
    @Test
    public void testBMPruneLikelihoodWithUpperMatrix10Species10Traits() {
        // tree
        treeStr = "(((t4:0.01700026351,t10:0.01700026351):0.9034244971,(t2:0.2519264674,t1:0.2519264674):0.6684982932):0.2901851714,((t9:0.4243929597,t5:0.4243929597):0.428956968,((t3:0.1329259849,t6:0.1329259849):0.5897285177,(t7:0.3444839115,t8:0.3444839115):0.3781705911):0.1306954251):0.3572600043):0.0;";
        spNames = "t4 t10 t2 t1 t9 t5 t3 t6 t7 t8";
        tree = new TreeParser(treeStr, false, false, true, 0);

        //trait values
        nTraits = 10;
        data = Arrays.asList(
                -2.28731878531972, 8.16016797286394, -1.60790363996585, -3.50375146318285, 2.4180310750393, 4.09620039236067, -5.17365781499908, 0.0830324263435673, -3.92697498145996, 1.23323556462526,
                -0.811378577256192, 8.18829611082131, -2.45648222243076, -4.4274634946743, 3.03356006715664, 3.3252839398938, -4.21928449405138, 0.150249836685309, -4.97807077472762, 1.03044047843148,
                13.6231428892787, -8.90573230225415, 4.21682244085758, 4.60332633627662, -0.716172061308149, 0.345250157267404, -0.749225468028548, -3.41701294175762, 0.994478698194049, -3.44063534931992,
                7.16152550409223, -12.1514214682735, 5.53944577979937, 7.14147083928846, -1.89982704153151, -2.16553733680629, -0.093001131322756, -3.32603366513231, 1.53951108454289, -2.64755502918416,
                -0.519423079492719, -5.77669193046342, -0.312497972990135, -1.59583025877346, 1.3446283692346, 2.99774158156863, 3.65399830020205, -1.69746045226351, -5.82006681351925, 0.84039105047492,
                -0.984130963068706, -3.0397054776387, 6.09565103724857, -0.840066990821478, 3.15088198649174, 1.77512813315887, 6.02316315859529, -2.73207211765912, -3.7025682293389, -0.571614944120547,
                6.92575877256499, -6.37840622885572, -1.76067360986415, -3.77388932553854, -0.453482143262971, 3.60403255588672, 1.3112428586143, -0.732024175817608, -3.42667016024356, 2.3074920977748,
                5.0957217584282, -5.25782891756779, 2.13664343352796, -2.64427540812803, -2.19363406695927, 2.25231513022396, 3.09693981602772, -1.64956110239375, -2.18351165272075, 2.14200618566428,
                5.26074139751262, 7.6536042206705, 5.13405841572984, -8.23652751442488, -0.431284906380698, 8.90155578270267, 0.0936214074138489, -3.00924775897672, -2.18342919122651, 1.66729459543154,
                4.83917711885396, 1.77623369001487, 3.99553867000827, -10.5516684942347, -3.257852292301, 6.47320619083183, 3.74483134602231, -3.62473323746156, -5.92016122137718, 3.14670114612779
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // BM model parameters
        sigmasq = new RealParameter(new Double[] {
                3.92444342735533, 3.02812661102414, 3.06556360045831, 2.81005713192123, 2.30077876608194, 4.64631271741818, 3.15615698487533, 1.50684016606037, 3.3548490285494, 2.35882255622108
        });
        correlation = new RealParameter(new Double[] {
                -0.714399955235422, -0.170907328370959, -0.172551347408444, -0.262309098150581, -0.695110504515469, -0.722387873101979, -0.533931801095605, -0.0680750994943082, -0.468054719269276, 0.715655430685729, -0.908337666653097, -0.11559985158965, 0.597849691286683, -0.75620148004964, 0.121895967517048, -0.58693722076714, -0.744936699513346, 0.506615728605539, 0.790090718306601, -0.251074448227882, 0.33023038925603, -0.810318678151816, -0.232060724403709, -0.451232710853219, 0.629280077759176, -0.102967317216098, 0.620128706097603, 0.624779019039124, 0.588684642221779, -0.120336624793708, 0.508950317278504, 0.258442263118923, 0.420364802703261, -0.998750453349203, -0.0493668518029153, -0.55976222967729, -0.24036692455411, 0.225542006548494, -0.296404181513935, -0.777729151304811, -0.512761054560542, 0.336111174896359, -0.164706440642476, 0.576391668058932, -0.794270711485296
        });
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation, "upperMatrix", true);

        // prune likelihood
        BMPruneLikelihood fpcm2 = new BMPruneLikelihood();
        fpcm2.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);
        double lnLk = fpcm2.calculateLogP();

        Assert.assertEquals(-205.515395500167, lnLk, EPSILON);
    }
}

