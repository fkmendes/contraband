package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.prunelikelihood.OUNodeMath;
import contraband.prunelikelihood.OUPruneLikelihood;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class OUPruneLikelihoodTest {
    final static double EPSILON = 1e-7;

    Tree tree;
    String treeStr;
    String spNames;
    int nTraits;
    List<Double> oneTraitValues;
    KeyRealParameter oneTraitData;

    OUNodeMath nodeMath;
    RealParameter rootValues;
    RealParameter sigma;
    RealParameter sigmae;
    RealParameter alpha;
    RealParameter theta;

    IntegerParameter colorAssignments;
    OUPruneLikelihood pcm;

    /*
     * (1) Ultrametric tree with 3 taxa, 2 continuous traits, 1 optima
     */
    @Test
    public void testOUPruneLkSmallTreeTwoTraitsOneOpt() {
        treeStr = "((A:3.0058179,B:3.0058179):4.350951,C:7.3567689);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "A B C";

        nTraits = 2;
        oneTraitValues = Arrays.asList(
                0.131394584822684, -0.19269144948362,
                -0.65643500381027,-0.373259036803319,
                -0.17399894787897, 1.27056761078824
        );
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {0.2, 1.3});
        sigma = new RealParameter(new Double[] {1.0, 0.6, 1.0});
        sigmae = new RealParameter(new Double[] {1.0, 0.3, 1.0});
        alpha = new RealParameter(new Double[] {2.0, 0.0, 0.0, 2.0});
        theta = new RealParameter(new Double[] {0.5, 0.5});

        nodeMath = new OUNodeMath();
        colorAssignments = new IntegerParameter();
        colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0 });
        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta, "sigma", sigma, "sigmae", sigmae, "root", rootValues,
                "optNr", 1, "optAssign", colorAssignments);

        pcm = new OUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(-7.58111239313721, logP, EPSILON);
    }

    /*
     * (2) Ultrametric tree with 4 taxa, 3 continuous traits, 1 optimum
     */
    @Test
    public void testOUPruneLkSmallTreeThreeTraitsOneOpt() {
        treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "sp1 sp2 sp3 sp4";

        nTraits = 3;
        oneTraitValues = Arrays.asList(
                2.90115170364898, 2.14268145872343, 3.2824264362529,
                2.62770319717086, 1.97031352144559, 4.28940778474205,
                3.96814434248716, 3.19423412470143, 3.49323574837925,
                8.88546166300816, 9.77485338099449, -3.65830173123419
        );
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        rootValues = new RealParameter(new Double[] {2.2, 1.3, 0.5});
        sigma = new RealParameter(new Double[] {0.35, 0.06, 0.06, 0.2, 0.4, 0.8});
        alpha = new RealParameter(new Double[] {1.8, -0.9, 1.2, -0.9, 1.0, 1.4, 1.2, 1.4, 1.6});
        theta = new RealParameter(new Double[] {5.0, 4.8, 1.2});
        nodeMath = new OUNodeMath();
        colorAssignments = new IntegerParameter();
        colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0});

        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta, "sigma", sigma, "sigmae", sigmae, "root", rootValues,
                "optNr", 1, "optAssign", colorAssignments);
        pcm = new OUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(-260.827821834006, logP, EPSILON);
    }

    /*
     * (3) Ultrametric tree with 4 taxa, 1 continuous trait, 3 optima
     * i.e. (1) testOUMVNLkOneTraitSmallTree3optCondOnRVEstimateRV in 'OUMVNLikelihoodOneTraitTest'
     *
     */
    @Test
    public void testOUPruneLkSmallTreeOneTraitThreeOpt(){
        treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "sp1 sp2 sp3 sp4";

        oneTraitValues = Arrays.asList(
                0.237649365136715,
                0.295018750722361,
                0.881225138279161,
                0.206222932069516);
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames);

        sigma = new RealParameter(new Double[]{ 0.006082604 });
        alpha = new RealParameter(new Double[]{ 7.390366 });
        rootValues = new RealParameter(new Double[]{ 3.182460e-10 });
        theta  = new RealParameter(new Double[]{ 0.206222932117995, 0.26633408087427, 0.88122539543514 });

        nodeMath = new OUNodeMath();
        colorAssignments = new IntegerParameter(new Integer[]{ 1, 1, 2, 0, 1, 0, 0 });
        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta, "sigma", sigma, "root", rootValues,
                "optNr", 3, "optAssign", colorAssignments, "upperMatrix", false);

        pcm = new OUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(9.9161059, logP, EPSILON);
    }

    /*
     * (4) Non-ultrametric tree with 17 taxa, 2 continuous traits, 1 optima
     */
    @Test
    public void testOUPruneLkNonUltraTreeTwoTraitsOneOpt(){
        treeStr = "(((((((t13:0.4153928797,t14:0.4153928797):0.08969920155,t12:0.5050920813):1.335489756,(t5:0.3553770905,t6:1.634213927):0.2063679102):0.4058592147,(t4:0.1414704627,(((t10:0.6854930613,t11:0.6854930613):0.5700674504,t8:1.255560512):0.3875149714,((t9:0.7462358178,(t15:0.1723162068,(t16:0.01956079107,t17:0.01956079107):0.1527554157):1.026534376):0.1565906523,t7:0.2508814445):0.2876342481):0.4626681167):0.1406974519):1.009880328,t3:2.630416914):0.5606009664,t2:0.0885847823):1.07569483,t1:0.7529914498);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "t13 t14 t12 t5 t6 t4 t10 t11 t8 t9 t15 t16 t17 t7 t3 t2 t1";

        nTraits  = 2;
        oneTraitValues = Arrays.asList(
                4.53371989048144, 3.78615589093471,
                4.50377056689979, 4.49781867425133,
                5.25059484540885, 4.58053549399816,
                4.88109957265796, 4.42012959766164,
                4.96705912624135, 4.67721690494978,
                5.70218989023257, 4.90098912630942,
                5.63375249441939, 5.01053920083377,
                4.88382283792164, 4.57084172030421,
                4.92673894105087, 4.17199941425637,
                5.05489624167441, 5.25741244815168,
                5.84698911383908, 5.60672464845884,
                5.84090830550805, 5.26950832346736,
                5.87871402744199, 5.4430403269941,
                5.43290478387275, 5.30423176486044,
                4.83073908334688, 4.63456351247378,
                6.16245952168041, 6.08003886940573,
                6.04360811109139, 6.18364476700973);
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        sigma = new RealParameter(new Double[]{ 0.299127897384494, 0.0402696667094893, 0.363185292789432 });
        alpha = new RealParameter(new Double[]{ 2.37692446163834, -1.6037269282734, -1.6037269282734, 2.00156092469912 });
        rootValues = new RealParameter(new Double[]{ 3.27435707319454, 10.4770161810953 });
        theta  = new RealParameter(new Double[]{ 4.9056017367667, 4.48321528780333 });

        nodeMath = new OUNodeMath();
        colorAssignments = new IntegerParameter(new Integer[]{ 0 });
        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta, "sigma", sigma, "root", rootValues,
                "optNr", 1, "optAssign", colorAssignments, "upperMatrix", false);

        pcm = new OUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(-6.04423062090828, logP, EPSILON);
    }

    /*
     * (5) Ultrametric tree with 4 taxa, 2 continuous traits, 2 optima
     */
    @Test
    public void testOUPruneLkSmallTreeTwoTraitsTwoOpt(){
        treeStr = "(((sp1:2,sp2:2):3,sp3:5):1,sp4:6);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "sp1 sp2 sp3 sp4";

        nTraits = 2;
        oneTraitValues = Arrays.asList(
                -0.852088501326528, -0.724973069661161, 0.0490989886775859, 0.478463675183718, 0.477307495634625, -0.0761208848808844, 0.826579895044279, 0.122281759527914
        );
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {0.2, 1.3});
        sigma = new RealParameter(new Double[] {1.0, 0.6, 1.0});
        alpha = new RealParameter(new Double[] {1.0, 1.2, 1.2, 1.0});
        theta = new RealParameter(new Double[] {4.2, 3.8, 1.0, 2.0});

        nodeMath = new OUNodeMath();
        colorAssignments = new IntegerParameter(new Integer[]{ 1, 1, 0, 0, 0, 0, 0});
        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta, "sigma", sigma, "root", rootValues,
                "optNr", 2, "optAssign", colorAssignments);

        pcm = new OUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(-95.4784626512006, logP, EPSILON);
    }

    /*
     * (6) Non-ultrametric tree with 5 taxa, 1 continuous trait, 1 optima
     *     Having 1 sampled ancestor
     */
    /*
    @Test
    public void testOUPruneLkSATreeThreeTraitsOneOpt(){
        treeStr = "((t1:1.0,t2:0.0):2.0,t3:3.0);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "t1 t2 t3";

        nTraits = 1;
        oneTraitValues = Arrays.asList(
                0.735789246422011, 0.793446047540599, 1.73244219768237
        );
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {0.908763274179301});
        sigma = new RealParameter(new Double[] {0.0551858804629736});
        alpha = new RealParameter(new Double[] {5.86615612063902e-12});
        theta = new RealParameter(new Double[] {18487533048.7273});

        nodeMath = new OUNodeMath();
        colorAssignments = new IntegerParameter(new Integer[]{0});
        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta, "sigma", sigma, "root", rootValues,
                "optNr", 1, "optAssign", colorAssignments, "upperMatrix", false);

        pcm = new OUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(-0.8070789215538116678772, logP, EPSILON);
    }
    */
}
