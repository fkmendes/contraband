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

public class DOULikelihoodTest {
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
    private RealParameter alpha;
    private RealParameter dAlpha;
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
    public void testDOULkSmallTreeTwoTraitsOneOpt() {
        treeStr = "((A:3.0,B:3.0):4.0,C:7.0);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "A B C";

        nTraits = 2;
        contTraitData = Arrays.asList(
                -8.35873276897194, -8.75604744861027,
                -10.4625464429157, -10.1448377936672,
                0.381165633761954, 2.85191670820294
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

        alpha = new RealParameter(new Double[] {-0.2, 0.1, 0.1, -0.2});
        dAlpha = new RealParameter(new Double[] {0.1, -0.5, -0.5, 0.1});
        theta = new RealParameter(new Double[] {0.5, 0.5});
        optNr = 1;
        optAssignment = new IntegerParameter(new Integer[] {0});
        modelParameter.initByName("alpha", alpha, "dAlpha", dAlpha, "theta", theta, "optNr", optNr, "optAssign", optAssignment, "trait", morphData);
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "popMatrix", sigmaEMatrix, "tree", tree, "rootValues", rootValues);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-15.7588578266153, logP, EPSILON);
    }

    @Test
    public void testDOUPruneLkLargeTreeThreeTraitsOneOpt(){
        treeStr = "(t1:0.1163718304,(((((t17:0.2129628541,t18:0.2129628541):0.07453562764,t15:0.2874984817):0.3977425382,t11:0.3512942379):1.036886306,((t10:0.1562518592,((t12:0.5782513076,((t27:0.02513471432,t28:0.02513471432):0.2426007684,t16:0.062600673):0.3105158249):0.03989647097,(t23:0.05547937392,t24:0.05547937392):0.5626684047):0.3678257258):0.2560823992,((t19:0.1913521399,t20:0.1913521399):0.04652416149,(t21:0.06584183368,t22:0.06584183368):0.1720344677):1.004179602):0.4800714219):0.7140543601,(t2:1.431522854,(((((t8:0.2399291188,t9:1.011775103):0.04211094753,t7:0.6530845034):0.4745764724,((t25:0.04090256846,t26:0.04090256846):0.8053102614,(t13:0.270469454,t14:0.2102897052):0.5508876351):0.9521176059):0.2546605777,(t4:1.43054298,(t5:0.5252044729,t6:0.1712156116):0.01835441935):0.6224480339):0.04785280739,t3:0.3209607885):0.1983279977):0.137009867):0.6474643778);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "t1 t17 t18 t15 t11 t10 t12 t27 t28 t16 t23 t24 t19 t20 t21 t22 t2 t8 t9 t7 t25 t26 t13 t14 t4 t5 t6 t3";

        nTraits = 3;
        contTraitData = Arrays.asList(
                0.0977528576448568, -0.66221556279696, 0.481559734811505, -5.61249083271616, -57.9536430080658, 92.7001918442323, -4.02777394670528, -58.0439749604025, 92.7259622054938, -4.33856896091605, -58.0655894467586, 92.8845315922523, -8.09933673447157, -35.5650551430264, 64.2351465701515, -0.62132256228082, -23.9893087263263, 38.2965777527808, -1.76599313724319, -72.2457286854258, 108.917614530889, -1.57194322861991, -72.8464951108445, 108.940044752068, -1.50425364275528, -72.949329585109, 108.99402392614, -0.326713005868558, -55.9045156270911, 83.7546234276814, -3.74189463653075, -71.896845577072, 109.280086937514, -3.38928781043747, -71.9263016940901, 109.231931382714, -2.16942764068026, -71.5986954636042, 107.85421526239, -0.537838986506029, -71.6066116336437, 107.900123979184, -2.03662303639931, -71.5550411689675, 108.199220633625, -1.75595973379703, -71.379467670017, 108.030519742416, 9.44679854988748, -32.8723336969521, 35.0016634124532, 0.4069609290985, -17.6263270500544, 28.3764826142553, -0.943323824280587, -50.7712046623882, 75.9896392344186, 1.36231437691018, -30.9490002223437, 44.7083811043788, -2.5697162075529, -82.8651744951717, 125.862785120815, -2.66740427488043, -82.9990168581494, 125.842864663444, 2.66007582765776, -84.5768309506808, 120.465731161128, 1.37843267812611, -78.3437173528689, 111.470955443813, 10.133155372353, -81.6201843644575, 102.917551312253, -2.47961951528567, -19.9830722627771, 35.1171043659274, 0.0572626873269174, -13.0305965741325, 22.0689750999815, 1.31402062235501, -6.39964882781452, 10.8105541109458
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {-0.378395144137122, -1.0222174833606, -0.709954899788967});
        sigmasq = new RealParameter(new Double[] {1.0, 0.4, 0.5});
        covariance = new RealParameter(new Double[] {0.2, 0.6, 0.25});
        sigmaMatrix.initByName("sigmasq", sigmasq, "covariance", covariance, "upperMatrix", true, "trait", morphData);
        alpha = new RealParameter(new Double[] {2.49371773910003, 1.06535764326987, 0.80644314985783, 1.48213018117132, 1.19980987349491, 1.7363067142808, 1.22750776736856, 3.91589870772859, 1.45027752473983});
        dAlpha = new RealParameter(new Double[] {-1.3536384305301, 1.31758089971532, 1.00093798084139, 2.90478784013022, 0.354806870270435, -0.783067421053588, 3.12420959893389, 1.88057079724842, 2.8317655401661});
        theta = new RealParameter(new Double[] {1.07598653072371, 1.84252850217709, -0.0936068635351809});
        optNr = 1;
        optAssignment = new IntegerParameter(new Integer[] {0});
        modelParameter.initByName("alpha", alpha, "dAlpha", dAlpha, "theta", theta, "optNr", optNr, "optAssign", optAssignment, "trait", morphData);
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "tree", tree, "rootValues", rootValues);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        OULikelihood likelihood = new OULikelihood();
        likelihood.initByName("params", modelParameter, "nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = likelihood.calculateLogP();
        Assert.assertEquals(-67.3184454033438, logP, EPSILON);
    }
}
