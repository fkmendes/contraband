package test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.prunelikelihood.MorphologicalData;
import contraband.prunelikelihood.ShrinkageParameter;
import contraband.utils.ShrinkageUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

public class ShrinkageParameterTest {
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> contTraitData;
    private final KeyRealParameter contTrait = new KeyRealParameter();
    private final MorphologicalData morphData = new MorphologicalData();
    private final ShrinkageParameter shrinkageParameter = new ShrinkageParameter();
    private int nSpecies;
    private double[] weight;
    private RealMatrix traitMat;
    private double delta;
    final static double EPSILON = 1e-7;

    @Test
    public void testShrinkageUtilsSixSpeciesFourTraits() {
        // initialize default weight vector
        nSpecies = 6;
        weight = new double [nSpecies];
        for (int i = 0; i < nSpecies; i++) {
            weight[i] = 1.0 / nSpecies;
        }

        // trait values
        traitMat = new Array2DRowRealMatrix(new double[][]{
                {1.43992069991891, -0.521710606400385, 1.10561637581434, -0.157612819068943},
                {0.348877438059058, 0.632563764185396, -0.157874694270553, -0.199347057660479},
                {-1.07819577781223, -0.466780487232386, -0.121903425587596, -0.735661003143472},
                { 1.08931843860522, -0.447764358092691, -0.200742989973232, -0.00327866708025785},
                {-1.14168814186063, -0.63123133899665, 0.01553830524637, -0.501053984275845},
                {-0.292510218946999, -1.56323463393269, -1.06283258902876, -0.75184508749932}
        });

        delta = ShrinkageUtils.estimateDelta(traitMat, weight);

        Assert.assertEquals(0.435026173002864, delta, EPSILON);
    }

    @Test
    public void testShrinkageUtilsFourSpeciesSixTraits(){
        // initialize default weight vector
        nSpecies = 4;
        weight = new double [nSpecies];
        for (int i = 0; i < nSpecies; i++) {
            weight[i] = 1.0 / nSpecies;
        }

        // trait values
        traitMat = new Array2DRowRealMatrix(new double[][]{
                {-0.420579569207178, 1.15032920533582, 0.0968521515822504, 1.8354220617026, 0.286678735881861, 0.836318736314219},
                {1.05096067831501, 0.158880815001583, 0.837912137759615, 1.47742100778513, -0.585856579077398, -0.366206661570247},
                {0.855288472381937, 0.524615741705375, 0.338739799123716, 0.0338961896720944, 1.39591803575264, -1.98905413542492},
                {0.745755760649241, -1.45262917873173, -2.46984489833647, 1.82339896004897, -0.900298166836283, 0.856360698349844}
        });

        delta = ShrinkageUtils.estimateDelta(traitMat, weight);

        Assert.assertEquals(0.676989791098336, delta, EPSILON);
    }

    @Test
    public void testShrinkageParameterSixSpeciesFourTraits(){
        // tree
        treeStr = "(t1:3.682105503,(t2:3.22892486,((t3:1.303684038,t4:1.303684038):1.056059251,(t5:1.773439676,t6:1.773439676):0.5863036128):0.8691815709):0.4531806428);";
        spNames = "t1 t2 t3 t4 t5 t6";
        tree = new TreeParser(treeStr, false, false, true, 0);
        nSpecies = 6;

        // trait values
        nTraits = 4;
        contTraitData =  Arrays.asList(
                1.43992069991891, -0.521710606400385, 1.10561637581434, -0.157612819068943,
                0.348877438059058, 0.632563764185396, -0.157874694270553, -0.199347057660479,
                -1.07819577781223, -0.466780487232386, -0.121903425587596, -0.735661003143472,
                1.08931843860522, -0.447764358092691, -0.200742989973232, -0.00327866708025785,
                -1.14168814186063, -0.63123133899665, 0.01553830524637, -0.501053984275845,
                -0.292510218946999, -1.56323463393269, -1.06283258902876, -0.75184508749932
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        shrinkageParameter.initByName("trait", morphData);

        delta = shrinkageParameter.getDelta();
        Assert.assertEquals(0.435026173002864, delta, EPSILON);
    }

    @Test
    public void testShrinkageParameterFourSpeciesSixTraits(){
        // tree
        treeStr = "(((t1:0.2724306467,t2:0.2724306467):0.1437224097,t3:0.4161530564):0.745661126,t4:1.161814182);";
        spNames = "t1 t2 t3 t4";
        tree = new TreeParser(treeStr, false, false, true, 0);
        nSpecies = 4;

        // trait values
        nTraits = 6;
        contTraitData =  Arrays.asList(
                -0.420579569207178, 1.15032920533582, 0.0968521515822504, 1.8354220617026, 0.286678735881861, 0.836318736314219,
                1.05096067831501, 0.158880815001583, 0.837912137759615, 1.47742100778513, -0.585856579077398, -0.366206661570247,
                0.855288472381937, 0.524615741705375, 0.338739799123716, 0.0338961896720944, 1.39591803575264, -1.98905413542492,
                0.745755760649241, -1.45262917873173, -2.46984489833647, 1.82339896004897, -0.900298166836283, 0.856360698349844
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        shrinkageParameter.initByName("trait", morphData);

        delta = shrinkageParameter.getDelta();
        Assert.assertEquals(0.676989791098336, delta, EPSILON);
    }

    @Test
    public void testShrinkageUtilsFiveSpeciesFiveTraits(){
        // initialize default weight vector
        nSpecies = 5;
        weight = new double [nSpecies];
        for (int i = 0; i < nSpecies; i++) {
            weight[i] = 1.0 / nSpecies;
        }

        // trait values
        traitMat = new Array2DRowRealMatrix(new double[][]{
                {-4.87477558734297, -4.59000477892991, 6.72451484591916, -5.18805142813354, -1.19197635519731},
                {5.39749477314946, -3.25048812813907, -2.07942745651652, -0.793602641268273, -3.57767696950484},
                {6.22112066579261, -3.07178215097587, 1.20017569174252, -0.945274137321795, -3.54976251898861},
                {9.92354205050096, -2.68844893108873, -4.07477870212002, 4.84454792985525, -4.23805923235145},
                {0.930842073808694, -3.47899525776201, -3.54978344328983, -0.959778634107881, -5.75213618335604}
        });

        delta = ShrinkageUtils.estimateDelta(traitMat, weight);

        Assert.assertEquals(0.528076427063438, delta, EPSILON);
    }

    @Test
    public void testShrinkageParameterFiveSpeciesFiveTraits(){
        // tree
        treeStr = "(t1:5.571467086,(t2:2.894831753,((t3:2.149684662,t4:2.149684662):0.5234919646,t5:2.673176626):0.2216551269):2.676635332);";
        spNames = "t1 t2 t3 t4 t5";
        tree = new TreeParser(treeStr, false, false, true, 0);
        nSpecies = 5;

        // trait values
        nTraits = 5;
        contTraitData =  Arrays.asList(
                -4.87477558734297, -4.59000477892991, 6.72451484591916, -5.18805142813354, -1.19197635519731,
                5.39749477314946, -3.25048812813907, -2.07942745651652, -0.793602641268273, -3.57767696950484,
                6.22112066579261, -3.07178215097587, 1.20017569174252, -0.945274137321795, -3.54976251898861,
                9.92354205050096, -2.68844893108873, -4.07477870212002, 4.84454792985525, -4.23805923235145,
                0.930842073808694, -3.47899525776201, -3.54978344328983, -0.959778634107881, -5.75213618335604
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        shrinkageParameter.initByName("trait", morphData);

        delta = shrinkageParameter.getDelta();
        Assert.assertEquals(0.528076427063438, delta, EPSILON);
    }

    @Test
    public void testShrinkageUtilsFiveSpeciesTwelveTraits() {
        // initialize default weight vector
        nSpecies = 5;
        weight = new double [nSpecies];
        for (int i = 0; i < nSpecies; i++) {
            weight[i] = 1.0 / nSpecies;
        }

        // trait values
        traitMat = new Array2DRowRealMatrix(new double[][]{
                {1.22757111274749, -4.03967634353308, 8.40538293948057, -0.470078247328622, -2.1369385249315, 3.32189880730566, 0.168065042257356, -5.77050415841039, -6.60432161327364, -4.65806967573443, 5.44154319382899, -1.17274756542552},
                {3.09296732034003, -3.34966577325527, 4.96474684786511, 3.76761092849908, -1.55222928028631, -3.91017429388793, 0.626780360664306, -6.01100918636374, -8.00397657760519, -5.7884403916855, 3.79613087411358, -1.59467300094332},
                {0.929144509718655, -9.77898860362989, 4.80306674778778, 0.754513402637935, -7.26307231302151, 1.26016300316014, -1.9963514296546, -0.315723773151065, -0.708401850965401, 1.53851144803393, 4.29442223556511, -2.01999240393342},
                {2.44255620923704, -8.99441090560193, 3.67037736568392, -0.0391320914929257, -7.84772464239177, 0.515919330618107, -4.9495368940019, 4.13781080682126, 0.406971964279643, -0.0754993151266581, 4.4985121251248, -2.80031955379674},
                {4.52472352123109, -7.92654085910766, 4.08090874957429, -0.0667941834152315, -6.37884031149302, 3.74704523898063, -4.74536527498033, 4.99883452612773, -7.32738746446273, 3.34218616234527, 5.54698788153555, -1.57790929469913}
        });

        delta = ShrinkageUtils.estimateDelta(traitMat, weight);

        Assert.assertEquals(0.543770621867578, delta, EPSILON);
    }

    @Test
    public void testShrinkageParameterFiveSpeciesTwelveTraits(){
        // tree
        treeStr = "((t1:2.743388858,t2:2.743388858):8.301125085,((t3:1.254641003,t4:1.254641003):1.691907801,t5:2.946548804):8.097965139);";
        spNames = "t1 t2 t3 t4 t5";
        tree = new TreeParser(treeStr, false, false, true, 0);
        nSpecies = 5;

        // trait values
        nTraits = 12;
        contTraitData =  Arrays.asList(
                1.22757111274749, -4.03967634353308, 8.40538293948057, -0.470078247328622, -2.1369385249315, 3.32189880730566, 0.168065042257356, -5.77050415841039, -6.60432161327364, -4.65806967573443, 5.44154319382899, -1.17274756542552, 3.09296732034003, -3.34966577325527, 4.96474684786511, 3.76761092849908, -1.55222928028631, -3.91017429388793, 0.626780360664306, -6.01100918636374, -8.00397657760519, -5.7884403916855, 3.79613087411358, -1.59467300094332, 0.929144509718655, -9.77898860362989, 4.80306674778778, 0.754513402637935, -7.26307231302151, 1.26016300316014, -1.9963514296546, -0.315723773151065, -0.708401850965401, 1.53851144803393, 4.29442223556511, -2.01999240393342, 2.44255620923704, -8.99441090560193, 3.67037736568392, -0.0391320914929257, -7.84772464239177, 0.515919330618107, -4.9495368940019, 4.13781080682126, 0.406971964279643, -0.0754993151266581, 4.4985121251248, -2.80031955379674, 4.52472352123109, -7.92654085910766, 4.08090874957429, -0.0667941834152315, -6.37884031149302, 3.74704523898063, -4.74536527498033, 4.99883452612773, -7.32738746446273, 3.34218616234527, 5.54698788153555, -1.57790929469913
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        shrinkageParameter.initByName("trait", morphData);

        delta = shrinkageParameter.getDelta();
        Assert.assertEquals(0.543770621867578, delta, EPSILON);
    }

    @Test
    public void testShrinkageUtilsTwelveSpeciesFiveTraits(){
        // initialize default weight vector
        nSpecies = 12;
        weight = new double [nSpecies];
        for (int i = 0; i < nSpecies; i++) {
            weight[i] = 1.0 / nSpecies;
        }

        // trait values
        traitMat = new Array2DRowRealMatrix(new double[][]{
                {-2.1897883394613, 0.0313359815175198, 5.02365697483978, -0.618578989331967, -3.01387472226807},
                {0.789127052367002, 0.479217272549513, 3.54132168390706, -2.09421282752067, -3.03657303810015},
                {2.80202830006273, 0.886130625647486, 4.59991939339537, -1.79290194413071, -4.81067400922374},
                {5.13578496477532, 1.22864060467125, 2.2813039275738, -1.03981260193573, -4.79665838428698},
                {-4.08812111642109, 0.240434594608957, -2.23929355990271, -2.28738549766977, -3.45188949026898},
                {-6.74180455928669, 3.96619049993955, 0.904152426360484, -3.35629321946682, 4.61310111815135},
                {-5.6100221124425, 3.30716488699528, 1.87786025723574, -2.27371691456047, 4.87854372712174},
                {1.04990602288301, 3.01608877474802, 1.24556137515722, -3.01784446762082, -1.54582728096936},
                {-3.44295019203676, 1.46502247056061, -2.52848949852231, 6.10229184573588, 3.54875058728228},
                {-5.71114447782265, 0.928413067163238, -2.58377857435135, 6.06015533050127, 1.63490091385978},
                {-0.906649100581506, 1.51753342890478, 0.170292544669327, 4.28370622379737, -1.89626789050161},
                {2.54048840314469, -1.84473872547521, 1.25933181381947, 2.46266563747088, -1.3022113394064}
        });

        delta = ShrinkageUtils.estimateDelta(traitMat, weight);

        Assert.assertEquals(0.330907701854751, delta, EPSILON);
    }

    @Test
    public void testShrinkageParameterTwelveSpeciesFiveTraits(){
        // tree
        treeStr = "((((t1:0.8630412891,t2:0.8630412891):1.932892233,(t3:1.071331611,t4:1.071331611):1.72460191):6.687482641,t5:9.483416163):0.7394673966,((t6:0.7348618051,t7:0.7348618051):6.796716618,((t8:6.73688377,(t9:2.987989272,t10:2.987989272):3.748894498):0.3196768994,(t11:2.102955104,t12:2.102955104):4.953605565):0.4750177538):2.691305136);";
        spNames = "t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12";
        tree = new TreeParser(treeStr, false, false, true, 0);
        nSpecies = 12;

        // trait values
        nTraits = 5;
        contTraitData =  Arrays.asList(
                -2.1897883394613, 0.0313359815175198, 5.02365697483978, -0.618578989331967, -3.01387472226807, 0.789127052367002, 0.479217272549513, 3.54132168390706, -2.09421282752067, -3.03657303810015, 2.80202830006273, 0.886130625647486, 4.59991939339537, -1.79290194413071, -4.81067400922374, 5.13578496477532, 1.22864060467125, 2.2813039275738, -1.03981260193573, -4.79665838428698, -4.08812111642109, 0.240434594608957, -2.23929355990271, -2.28738549766977, -3.45188949026898, -6.74180455928669, 3.96619049993955, 0.904152426360484, -3.35629321946682, 4.61310111815135, -5.6100221124425, 3.30716488699528, 1.87786025723574, -2.27371691456047, 4.87854372712174, 1.04990602288301, 3.01608877474802, 1.24556137515722, -3.01784446762082, -1.54582728096936, -3.44295019203676, 1.46502247056061, -2.52848949852231, 6.10229184573588, 3.54875058728228, -5.71114447782265, 0.928413067163238, -2.58377857435135, 6.06015533050127, 1.63490091385978, -0.906649100581506, 1.51753342890478, 0.170292544669327, 4.28370622379737, -1.89626789050161, 2.54048840314469, -1.84473872547521, 1.25933181381947, 2.46266563747088, -1.3022113394064
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        shrinkageParameter.initByName("trait", morphData);

        delta = shrinkageParameter.getDelta();
        Assert.assertEquals(0.330907701854751, delta, EPSILON);
    }

    @Test
    public void testShrinkageParameterWeightsOnTenSpeciesSixTraits(){
        // tree
        treeStr = "((((t1:5.001472856,(t9:0.8671472739,t10:0.8671472739):4.134325582):1.169748763,(t2:4.572865663,t3:4.572865663):1.598355956):1.028509947,(t6:4.061841536,(t7:1.217047531,t8:1.217047531):2.844794005):3.137890029):4.120313239,(t4:4.190072813,t5:4.190072813):7.129971992);";
        spNames = "t1 t9 t10 t2 t3 t6 t7 t8 t4 t5";
        tree = new TreeParser(treeStr, false, false, true, 0);
        nSpecies = 10;

        // trait values
        nTraits = 6;
        contTraitData =  Arrays.asList(
                -6.77771511714266, 5.6412376528486, 16.7787303972049, -0.808995434701036, 8.84439473245568, -6.03304843423229, 0.0240147807850825, -3.60199374802833, 6.70746847906576, 2.77219554352954, 5.48043562292939, -1.63861222545275, 1.21029257061537, -2.41207663256475, 5.21860888870714, 4.77394782636081, 8.22273655620462, -2.37967810465363, 4.45362620311193, -1.09471538833579, 0.259051733358624, 3.97746835296685, -2.41635825319525, -0.412562158053475, -4.28311743997789, -8.96041442159479, 11.5338063671715, 4.43082991266315, -0.395911017960692, -4.00877551895379, 2.91158824852438, -1.20271538091808, 11.2447167807583, 3.37338070586235, 10.2052105388768, -6.37663006990716, 5.49823972598282, -0.968407499325085, 6.33107250956359, 2.69726849358064, 2.40651903200087, -4.31901386089944, -1.3847760541425, 8.43167800720544, 11.6639685126154, -1.76960791902215, 4.59559070589586, -5.99933056361048, 7.69273157577489, 4.67428183396886, 3.37204425624012, 1.3060274229292, 2.07136812499068, -1.40548918534553, 2.68306230316882, 15.8942564699219, -5.05827986921687, 2.9907735251968, 3.96583930851246, 6.92390734541384
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);

        Double[] weight = new Double []{0.1, 0.1, 0.2, 0.05, 0.20, 0.05, 0.05, 0.15, 0.05, 0.05};

        RealParameter traitWeight = new RealParameter(weight);
        shrinkageParameter.initByName("trait", morphData, "weight", traitWeight);

        delta = shrinkageParameter.getDelta();
        Assert.assertEquals(0.90016696742332, delta, EPSILON);
    }

}
