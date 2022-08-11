package contraband.test;

import beast.base.evolution.tree.Tree;
import org.junit.Assert;
import org.junit.Test;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeParser;
import contraband.mvnlikelihood.BMMVNShiftLikelihoodOneTrait;
import contraband.clock.RateCategoryClockModel;
import contraband.clock.TreeToVCVMat;
// import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

/**
 * @author Fabio K. Mendes
 */

/*
 * This class contains unit tests for BMMVNShiftLikelihoodOneTrait, which allows for rate shifts (without and with coalescent correction)
 */
public class BMMVNShiftLikelihoodOneTraitTest {

    final static double EPSILON = 1e-6;

    Tree myTree;
    String treeStr;
    String spNames;
    RealParameter colorValues;
    IntegerParameter colorAssignments;
    RateCategoryClockModel rcc;
    TreeToVCVMat colors;
    Double[] rootValueVectorInput;
    RealParameter rootValue;
    List<Double> oneTraitValues;
    RealParameter oneTraitData;
    // KeyRealParameter oneTraitData;
    BMMVNShiftLikelihoodOneTrait bmLk;

    /*
     * (1) Small tree, one trait, one rate, without and with root edge (of 1.0). Should match/reduce to BMMVNLikelihoodOneTraitTest.
     */
    @Test
    public void testBMMVNShiftLkOneTraitSmallTree() {
        // tree
        treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
        myTree = new TreeParser(treeStr, false, false, true, 0);

        // VCV Mat
        colorValues = new RealParameter(new Double[]{ 0.2704762 });
        colorAssignments = new IntegerParameter(new Integer[]{0, 0, 0, 0, 0});
        // IntegerParameter rootEdgeColorAssignment = new IntegerParameter(new Integer[] { 0 });
        rcc = new RateCategoryClockModel();
        rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

        colors = new TreeToVCVMat();
        colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);
        TreeToVCVMat colors2 = new TreeToVCVMat();
        colors2.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", true, "rootEdgeLength", 1.0);

        /* DEPRECATED (original implementation of a sort of morph clock model) */
        // ColorManager colors = new ColorManager();
        // colors.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
        // ColorManager colors2 = new ColorManager();
        // colors2.initByName("nTraits", 1, "nColors", 1, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", true, "rootEdgeLength", 1.0, "rootEdgeColorAssignment", rootEdgeColorAssignment);

        // initializing data
        spNames = "sp1 sp2 sp3";
        oneTraitValues = Arrays.asList(4.1, 4.5, 5.9);
        oneTraitData = new RealParameter();
        // oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

        // root value vector
        rootValueVectorInput = new Double[] { 4.985714 };
        rootValue = new RealParameter(rootValueVectorInput);

        // likelihood
        bmLk = new BMMVNShiftLikelihoodOneTrait();
        bmLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
        double lnLk = bmLk.calculateLogP(); // no root edge

        bmLk = new BMMVNShiftLikelihoodOneTrait();
        bmLk.initByName("tree", myTree, "rateManager", colors2, "rootValue", rootValue, "oneTraitData", oneTraitData);
        double lnLk2 = bmLk.calculateLogP();

        Assert.assertEquals(-3.191339, lnLk, EPSILON);
        Assert.assertEquals(-3.577933, lnLk2, EPSILON);
    }

    /*
     * (2) Large ultrametric tree, one trait, one rate, no shifts.
     */
    @Test
    public void testBMMVNShiftLkOneTraitLargeTree() {
        // tree
        treeStr = "((((((t40:5.88018515,t30:5.88018515):31.84236817,(((t43:2.534803909,t25:2.534803909):12.53459081,t1:15.06939472):15.93712814,t4:31.00652286):6.716030459):4.293998055,t7:42.01655137):47.42017615,((t34:39.94317937,t13:39.94317937):42.87554852,t39:82.81872788):6.617999642):2.795514817,(((t19:31.28909285,t50:31.28909285):21.85222958,t9:53.14132243):28.5184562,t26:81.65977864):10.57246371):7.767757656,((t24:69.64568591,((((t5:7.026526759,t27:7.026526759):10.01835474,t45:17.0448815):30.78336706,t8:47.82824856):20.30558317,(((t41:5.397952755,t37:5.397952755):27.52253045,t42:32.9204832):24.003255,t11:56.9237382):11.21009353):1.511854173):12.42349092,(((t16:69.96444119,((t32:24.5178924,t14:24.5178924):37.67637018,(t2:16.45170245,t49:16.45170245):45.74256014):7.7701786):0.3095854321,((t33:30.05535248,((t23:13.30542076,t6:13.30542076):9.449794386,((t48:12.40052663,t21:12.40052663):2.499221124,t44:14.89974775):7.855467399):7.30013733):37.46376453,((t3:46.09553051,t29:46.09553051):5.973201324,t20:52.06873183):15.45038518):2.75490961):6.151101537,((((t28:24.94062793,t46:24.94062793):5.652833062,(t22:27.99970449,(t35:7.370771106,t38:7.370771106):20.62893338):2.593756503):10.27347673,(t17:25.27026498,t10:25.27026498):15.59667274):21.4744838,((t47:6.373725141,t36:6.373725141):45.36629128,(((t15:0.03406798076,t31:0.03406798076):2.516336885,t18:2.550404866):27.30416891,t12:29.85457377):21.88544265):10.60140509):14.08370664):5.644048672):17.93082317):0;";
        myTree = new TreeParser(treeStr, false, false, true, 0);

        // VCV Mat
        colorValues = new RealParameter(new Double[] { 0.06650558 });
        colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
        rcc = new RateCategoryClockModel();
        rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

        colors = new TreeToVCVMat();
        colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);

        // initializing data
        oneTraitValues = Arrays.asList(-1.98102089083783, 0.425095619275142, -2.64278670721857, 4.67301038174566, 0.980837136224631, -1.50747992926933, 0.291651433101693, 1.38909936714432, 4.03813926115063, 2.67745449448204, 6.39102061816721, 6.40419053691823, 2.20955411315545, 0.772842105471847, 3.55734391167665, 0.498071400234979, 0.550315407342643, 3.84375731077216, 3.18235300029237, 3.91668458860299, 5.2164188658234, 3.39905904957137, 2.69218348604267, -3.93463904654453, 2.03170767420157, 0.996489601169867, 0.742277069809021, -1.83602057962268, -1.01090704411752, 3.03382850811265, 5.43622150183381, 1.74925732548084, 3.5197105768139, 3.19459541434014, -0.457078134951209, 0.193865716678944, -0.088604135083334, 0.12961684230114, 1.65864752553039, 1.73675049380547, 1.25491657987277, 1.37659283253291, 2.25570866291172, 1.14726492498239, 0.25210113076106, 0.776384238776296, 1.31397115366671, -0.0419698505695754, -0.21193756283334, -0.138788530248324);
        spNames = "t39 t26 t9 t7 t34 t13 t19 t50 t4 t1 t40 t30 t43 t25 t16 t24 t11 t20 t8 t3 t29 t42 t33 t12 t22 t17 t10 t28 t46 t32 t14 t45 t2 t49 t44 t23 t6 t48 t21 t35 t38 t5 t27 t47 t36 t41 t37 t18 t15 t31";
        oneTraitData = new RealParameter();
        // oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

        // root value vector
        rootValueVectorInput = new Double[] { 1.212958 };
        rootValue = new RealParameter(rootValueVectorInput);

        // likelihood
        BMMVNShiftLikelihoodOneTrait bmLk = new BMMVNShiftLikelihoodOneTrait();
        bmLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
        double lnLk = bmLk.calculateLogP();

        Assert.assertEquals(-94.47619, lnLk, 1e-5);
    }

    /*
     * (3) Large NON-ultrametric tree, one trait, one rate, no shifts.
     */
    @Test
    public void testBMMVNShiftLkOneTraitLargeTreeNonUltra() {
        // tree
        treeStr = "(((((t35:0.1,t32:0.1):0.1,t10:0.1):0.1,t18:0.1):0.1,(((t47:0.1,t9:0.1):0.1,(t43:0.1,t38:0.1):0.1):0.1,(((((t20:0.1,t14:0.1):0.1,t19:0.1):0.1,(t24:0.1,(((t50:0.1,t8:0.1):0.1,t25:0.1):0.1,(t12:0.1,t5:0.1):0.1):0.1):0.1):0.1,t37:0.1):0.1,(t42:0.1,(t13:0.1,t41:0.1):0.1):0.1):0.1):0.1):0.1,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0;";
        myTree = new TreeParser(treeStr, false, false, true, 0);

        // VCV Mat
        colorValues = new RealParameter(new Double[]{0.0676543});
        colorAssignments = new IntegerParameter(new Integer[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
        rcc = new RateCategoryClockModel();
        rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

        colors = new TreeToVCVMat();
        colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);

        // initializing data
        oneTraitValues = Arrays.asList(-0.140079330558364, -0.0874251270522149, -0.194010550881256, 0.129301918296212, 0.0774210684049538, -0.047482768983324, 0.0956948735831467, 0.162500177959304, 0.18863127851286, 0.125036602967829, 0.350150590083319, 0.272741469616044, 0.116302901473441, -0.150985676720302, 0.236398643978856, 0.18383828082074, 0.184437478591962, 0.306155108012576, 0.308705544566368, 0.356922129276011, 0.437517619840372, 2.47333215098386, 0.70777639704284, -2.62662966633079, 0.454522834320791, -0.132242095134662, -0.439104186981486, -3.92360496759795, -2.41829387918679, 0.101516116711762, 2.6441683780108, 2.13821582254339, 2.54220446636924, -0.520807491795214, -4.29749916896601, -2.68665542778722, -3.01604647309455, -2.57037860992468, 2.0262053690181, 2.69837639068034, 1.01248730643342, 0.662522088866049, 1.70637323736846, 2.25896738991068, 2.20128394055548, -1.06070652973863, 1.20670535726105, 2.35323217229207, 1.91628626588915, 2.92725041007494);
        spNames = "t35 t32 t10 t18 t47 t9 t43 t38 t20 t14 t19 t24 t50 t8 t25 t12 t5 t37 t42 t13 t41 t34 t4 t36 t7 t29 t22 t46 t40 t28 t33 t21 t26 t48 t39 t2 t44 t23 t11 t49 t45 t31 t16 t30 t1 t17 t15 t3 t27 t6";
        oneTraitData = new RealParameter();
        // oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

        // root value vector
        rootValueVectorInput = new Double[]{0.1007117};
        rootValue = new RealParameter(rootValueVectorInput);

        // likelihood
        BMMVNShiftLikelihoodOneTrait bmLk = new BMMVNShiftLikelihoodOneTrait();
        bmLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
        double lnLk = bmLk.calculateLogP();

        Assert.assertEquals(-42.22841, lnLk, 1e-5);
    }

    /*
     * (4) Small ultrametric tree, one trait, two rates (one shift).
     */
    @Test
    public void testBMMVNShiftLkOneTraitSmallTreeTwoRates() {
        // tree
        treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
        myTree = new TreeParser(treeStr, false, false, true, 0);

        // VCV Mat
        colorValues = new RealParameter(new Double[] { 0.05057867, 3.360241 });
        colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 1, 1, 1 });
        rcc = new RateCategoryClockModel();
        rcc.initByName("nCat", 2, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

        colors = new TreeToVCVMat();
        colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);

        // initializing data
        oneTraitValues = Arrays.asList(-2.53718502574816, -2.85562629168723, 1.79661600241838);
        spNames = "sp1 sp2 sp3";
        oneTraitData = new RealParameter();
        // oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

        // root value vector
        rootValueVectorInput = new Double[] { -1.191236 };
        rootValue = new RealParameter(rootValueVectorInput);

        // likelihood
        bmLk = new BMMVNShiftLikelihoodOneTrait();
        bmLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
        double lnLk = bmLk.calculateLogP();

        Assert.assertEquals(-4.673609, lnLk, EPSILON);
    }

    /*
     * (5) Large ultrametric tree, one trait, three rates (two shifts).
     */
    @Test
    public void testBMMVNShiftLkOneTraitLargeTreeThreeRates() {
        // tree
        treeStr = "(((((t35:2.336518061,t32:2.336518061):28.95257479,t10:31.28909285):8.654086516,t18:39.94317937):52.28906298,(((t47:31.00652286,t9:31.00652286):50.20634817,(t43:15.06939472,t38:15.06939472):66.14347631):10.61662549,(((((t20:20.94406932,t14:20.94406932):28.09437292,t19:49.03844224):31.16991698,(t24:54.88723469,(((t50:2.534803909,t8:2.534803909):35.18774941,t25:37.72255332):15.41876911,(t12:42.01655137,t5:42.01655137):11.12477106):1.745912255):25.32112453):2.610368667,t37:82.81872788):6.617999642,(t42:81.65977864,(t13:5.88018515,t41:5.88018515):75.77959349):7.776948892):2.392768999):0.4027458181):7.767757656,((t34:80.73867518,((t4:14.89974775,t36:14.89974775):7.855467399,t7:22.75521515):57.98346003):16.48666894,((((((((t29:32.9204832,t22:32.9204832):13.17504731,t46:46.09553051):1.732718052,t40:47.82824856):14.51317295,(t28:29.85457377,((t33:6.373725141,t21:6.373725141):1.191235246,t26:7.564960387):22.28961339):32.48684774):5.177695495,t48:67.51911701):2.445324178,(t39:56.9237382,((t2:5.876590264,t44:5.876590264):19.06403767,t23:24.94062793):31.98311027):13.04070299):0.3095854321,(((t11:13.30542076,t49:13.30542076):14.69428372,t45:27.99970449):1.437902517,t31:29.43760701):40.83641961):11.48412211,((((t16:30.59346099,(t30:0.03406798076,t1:0.03406798076):30.55939301):21.47527084,(t17:50.41024027,t15:50.41024027):1.658491566):14.63237622,(t3:10.35007739,t27:10.35007739):56.35103066):2.944577857,t6:69.64568591):12.11246283):15.46719539):2.774655878):0.0;";
        myTree = new TreeParser(treeStr, false, false, true, 0);

        // VCV Mat
        colorValues = new RealParameter(new Double[] { 0.006850499, 0.08316839, 0.2736696 });

        // first 5 rows are tips
        colorAssignments = new IntegerParameter(new Integer[] {
                0, 0, 0, 2, 2, 2, 0, 0, 0, 0,
                2, 1, 2, 1, 1, 1, 2, 2, 1, 0,
                1, 1, 0, 0, 0, 0, 1, 0, 0, 0,
                2, 0, 1, 0, 1, 2, 2, 0, 1, 0,
                1, 0, 1, 0, 2, 2, 0, 0, 2, 0,

                0, 0, 0, 0, 0, 0, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 0, 0,
                0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0
        });
        rcc = new RateCategoryClockModel();
        rcc.initByName("nCat", 3, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

        colors = new TreeToVCVMat();
        colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);
        // ColorManager colors = new ColorManager();
        // colors.initByName("tree", myTree, "nTraits", 1, "nColors", 3, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);

        // initializing data
        oneTraitValues = Arrays.asList(-4.24865266388791, 1.12694466883781, -6.02979129548445, 9.62359673885322, 0.17990827089315, -5.46776437171474, 0.486983190714245, 2.08585804088066, 0.753683973121965, -0.242019913226091, 0.946577159488485, 6.22141276151274, 2.91113662785022, -2.0374126462967, -1.10293312619706, -0.179708562978063, -0.198701248227853, 4.58084601841762, 5.64277550994041, 0.919256884182315, 1.09524719371436, 0.782136300731067, 0.0962378659531526, -5.33631896191408, 1.4905959236315, -0.0744248296375519, -0.177479586475342, -4.26961572362966, -2.86137408169925, -0.974994208280936, 1.41858116564816, -0.187048179444733, 0.282668872536788, -0.244597214204735, -1.19517123879371, -1.10037443643529, -0.261483237873604, -0.207560694831451, 0.348278374786552, -0.168277479250151, -0.406155392866917, -0.3594972601634, -0.0234868583801822, 1.00637046502527, -0.000469931947963076, -0.329563226090097, 0.043226745691999, 0.739816593725853, -0.238002140615943, -0.214951085346446);
        spNames = "t37 t42 t24 t19 t12 t5 t18 t25 t10 t47 t9 t20 t14 t43 t38 t13 t41 t50 t8 t35 t32 t34 t6 t48 t39 t17 t15 t40 t46 t29 t22 t16 t28 t31 t45 t23 t7 t4 t36 t11 t49 t3 t27 t26 t33 t21 t2 t44 t30 t1";
        oneTraitData = new RealParameter();
        // oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

        // root value vector
        rootValueVectorInput = new Double[] { -0.0501181 };
        rootValue = new RealParameter(rootValueVectorInput);

        // likelihood
        bmLk = new BMMVNShiftLikelihoodOneTrait();
        bmLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
        double lnLk = bmLk.calculateLogP();

        Assert.assertEquals(-80.18669, lnLk, 1e-5);
    }

    /*
     * (6) Large NON-ultrametric tree with sampled ancestors, one trait, one rate.
     */
    @Test
    public void testBMMVNShiftLkOneTraitLargeTreeNonUltraSampledAnc() {
        // tree
        treeStr = "((((t16_1:53.23751153,t5_1:53.23751153):8.655770254,((t2_3:8.807060933,t2_2:0):35.79300727,t2_1:0):17.29321358):42.65570341,((((t13_4:37.22097358,t13_3:0):16.65221927,t13_2:0):5.742661592,t13_1:0):32.39904669,((t7_1:29.60227551,t8_1:29.60227551):52.89177712,t28_1:0):9.520848496):12.53408407):67.09599568,((((t15_1:14.38089478,t3_1:0.5275386099):115.9798433,((t1_1:29.18219406,(((t10_3:2.464982927,t10_2:0):1.709987849,t10_1:0):22.8996627,t14_1:2.142260068):42.13119303):28.80221822,t11_1:8.156125409):32.35269336):29.28086375,t19_1:0):7.958804008,((((t4_1:110.1538492,((t6_1:10.73085343,t12_1:16.35390184):26.57130888,t9_1:29.24728433):67.22863851):17.9062182,t20_3:0):5.616707849,t20_2:0):1.720902578,t20_1:0):32.20272799):4.044575033):20.90108373;";
        myTree = new TreeParser(treeStr, false, false, true, 0);

        // VCV Mat
        colorValues = new RealParameter(new Double[] { 0.06715078 });
        colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
        rcc = new RateCategoryClockModel();
        rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

        colors = new TreeToVCVMat();
        colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);

        // initializing data
        oneTraitValues = Arrays.asList(-5.05744071356103, -1.11052196957346, -0.575427737395457, -0.202953693279079, 1.14969399531796, -1.64202406456342, 1.79417473949251, 0.073438714146757, 0.726726641446067, 3.93444954633643, -0.491682932895839, -1.74464882096195, -0.0780925842656348, 4.81261164374873, 4.49854934867917, 3.04106820281915, -0.803750712288891, 0.523401664662935, -0.456279043339767, 1.2741332640877, 2.87956906600202, 3.91356152320114, 3.59597823070803, -0.389239530235073, -0.63498836297618, -0.789764345711418, -1.20679458892863, -1.48400428008497);
        spNames = "t4_1 t12_1 t2_3 t7_1 t13_4 t15_1 t16_1 t5_1 t8_1 t10_3 t6_1 t9_1 t3_1 t14_1 t1_1 t11_1 t2_1 t2_2 t13_1 t13_2 t13_3 t10_1 t10_2 t19_1 t20_1 t20_2 t20_3 t28_1";
        oneTraitData = new RealParameter();
        // oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

        // root value vector
        rootValueVectorInput = new Double[] { -0.4990207 };
        rootValue = new RealParameter(rootValueVectorInput);

        // likelihood
        bmLk = new BMMVNShiftLikelihoodOneTrait();
        bmLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
        double lnLk = bmLk.calculateLogP();

        Assert.assertEquals(-48.54424, lnLk, 1e-5);
    }

    /*
     * (7) Large NON-ultrametric tree with sampled ancestors, no fossils, one trait, one rate.
     */
    @Test
    public void testBMMVNShiftLkOneTraitLargeTreeNonUltraSampledAncNoFossils() {
        // tree
        treeStr = "((((t16_1:53.23751153,t5_1:53.23751153):8.655770254,((t2_3:8.807060933,t2_2:0):35.79300727,t2_1:0):17.29321358):42.65570341,((((t13_4:37.22097358,t13_3:0):16.65221927,t13_2:0):5.742661592,t13_1:0):32.39904669,((t7_1:29.60227551,t8_1:29.60227551):52.89177712,t28_1:0):9.520848496):12.53408407):67.09599568,(((t15_1:130.3607381,((t10_3:2.464982927,t10_2:0):1.709987849,t10_1:0):126.1857673):29.28086375,t19_1:0):7.958804008,((((t4_1:110.1538492,t12_1:110.1538492):17.9062182,t20_3:0):5.616707849,t20_2:0):1.720902578,t20_1:0):32.20272799):4.044575033):20.90108373;";
        myTree = new TreeParser(treeStr, false, false, true, 0);

        // VCV Mat
        colorValues = new RealParameter(new Double[] { 0.07044749 });
        colorAssignments = new IntegerParameter(new Integer[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
        rcc = new RateCategoryClockModel();
        rcc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", myTree);

        colors = new TreeToVCVMat();
        colors.initByName("branchRateModel", rcc, "tree", myTree, "coalCorrection", false);

        // initializing data
        oneTraitValues = Arrays.asList(-5.05744071356103, -1.11052196957346, -0.575427737395457, -0.202953693279079, 1.14969399531796, -1.64202406456342, 1.79417473949251, 0.073438714146757, 0.726726641446067, 3.93444954633643, 0.750504711256288, -0.587613596378262, 2.75170560121156, 3.31308844229513, 3.23260100181807, 3.97279694209149, 3.75757747655993, 1.43752567783899, -1.53225408723805, -0.667885599149076, 0.266271594649991, 0.98828295113675);
        spNames = "t4_1 t12_1 t2_3 t7_1 t13_4 t15_1 t16_1 t5_1 t8_1 t10_3 t2_1 t2_2 t13_1 t13_2 t13_3 t10_1 t10_2 t19_1 t20_1 t20_2 t20_3 t28_1";
        oneTraitData = new RealParameter();
        // oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

        // root value vector
        rootValueVectorInput = new Double[] { 0.8925689 };
        rootValue = new RealParameter(rootValueVectorInput);

        // likelihood
        bmLk = new BMMVNShiftLikelihoodOneTrait();
        bmLk.initByName("tree", myTree, "rateManager", colors, "rootValue", rootValue, "oneTraitData", oneTraitData);
        double lnLk = bmLk.calculateLogP();

        Assert.assertEquals(-38.44236, lnLk, 1e-5);
    }
}
