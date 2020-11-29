package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import contraband.prunelikelihood.OUPruneLikelihood;
import org.apache.commons.math3.linear.*;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.prunelikelihood.OUNodeMath;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

public class OUNodeMathTest {
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
    OUPruneLikelihood pcm;
    IntegerParameter colorAssignments;

    @Test
    public void testParameterInSolvedExample() {
        treeStr = "((A:3.0058179,B:3.0058179):4.350951,C:7.3567689);";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "A B C";

        nTraits = 2;
        oneTraitValues = Arrays.asList(
                4.53371989048144, 4.39037816144705,
                4.64369804601356, 4.50915779102972,
                6.22349219457299, 6.60199054481913
        );
        oneTraitData = new KeyRealParameter();
        oneTraitData.initByName("value", oneTraitValues, "keys", spNames, "minordimension", nTraits);

        // OU model parameters
        rootValues = new RealParameter(new Double[] {0.2, 1.3});
        sigma = new RealParameter(new Double[] {1.36, 0.6, 1.0});
        sigmae = new RealParameter(new Double[] {1.0, 0.3, 1.0});
        alpha = new RealParameter(new Double[] {1.0, 3.0, 2.0, 4.0});
        theta = new RealParameter(new Double[] {0.5, 0.5});

        colorAssignments = new IntegerParameter();
        colorAssignments = new IntegerParameter(new Integer[]{ 0, 0, 0, 0, 0, 0 });

        nodeMath = new OUNodeMath();
        nodeMath.initByName("traits", oneTraitData, "alpha", alpha, "theta", theta,
                "sigma", sigma, "sigmae", sigmae,
                "root", rootValues,"optNr", 1, "optAssign", colorAssignments,
                "upperMatrix", false);

        pcm = new OUPruneLikelihood();
        pcm.initByName("tree", tree, "traits", oneTraitData, "nodeMath", nodeMath);

        double logP = pcm.calculateLogP();
        Assert.assertEquals(-56.447173900514, logP, EPSILON);

        // Matrix A
        RealMatrix a1RM = nodeMath.getAMatForNode(0);
        Double[] a1Mat = new Double[] {-0.10070390510264365, -0.1078030493354475, -0.32069931328358015};
        Double[] a1RMDouble = new Double[] {a1RM.getEntry(0,0), a1RM.getEntry(0,1), a1RM.getEntry(1,1)};
        Assert.assertArrayEquals(a1Mat, a1RMDouble);

        RealMatrix a2RM = nodeMath.getAMatForNode(1);
        Double[] a2Mat = new Double[] {-0.10070390510264365, -0.1078030493354475, -0.32069931328358015};
        Double[] a2RMDouble = new Double[] {a2RM.getEntry(0,0), a2RM.getEntry(0,1), a2RM.getEntry(1,1)};
        Assert.assertArrayEquals(a2Mat, a2RMDouble);

        RealMatrix a3RM = nodeMath.getAMatForNode(2);
        Double[] a3Mat = new Double[] {-0.06400003300242779, -0.1355806427516654, -0.29967715262733896};
        Double[] a3RMDouble = new Double[] {a3RM.getEntry(0,0), a3RM.getEntry(0,1), a3RM.getEntry(1,1)};
        Assert.assertArrayEquals(a3Mat, a3RMDouble);

        RealMatrix a4RM = nodeMath.getAMatForNode(3);
        Double[] a4Mat = new Double[] {-0.6230777350650233, -1.3268889621850428, -2.948141249636984};
        Double[] a4RMDouble = new Double[] {a4RM.getEntry(0,0), a4RM.getEntry(0,1), a4RM.getEntry(1,1)};
        Assert.assertArrayEquals(a4Mat, a4RMDouble);

        // Matrix C
        RealMatrix c1RM = nodeMath.getCMatForNode(0);
        Double[] c1Mat = new Double[] {-0.3757163513481878, 0.25779426672220895, -0.17688312929786404};
        Double[] c1RMDouble = new Double[] {c1RM.getEntry(0,0), c1RM.getEntry(0,1), c1RM.getEntry(1,1)};
        Assert.assertArrayEquals(c1Mat, c1RMDouble);

        RealMatrix c2RM = nodeMath.getCMatForNode(1);
        Double[] c2Mat = new Double[] {-0.3757163513481878, 0.25779426672220895, -0.17688312929786404};
        Double[] c2RMDouble = new Double[] {c2RM.getEntry(0,0), c2RM.getEntry(0,1), c2RM.getEntry(1,1)};
        Assert.assertArrayEquals(c2Mat, c2RMDouble);

        RealMatrix c3RM = nodeMath.getCMatForNode(2);
        Double[] c3Mat = new Double[] {-0.36977746373559073, 0.25371935362506826, -0.17408716516578368};
        Double[] c3RMDouble = new Double[] {c3RM.getEntry(0,0), c3RM.getEntry(0,1), c3RM.getEntry(1,1)};
        Assert.assertArrayEquals(c3Mat, c3RMDouble);

        RealMatrix c4RM = nodeMath.getCMatForNode(3);
        Double[] c4Mat = new Double[] {-0.38493740617556593, 0.2641212065640522, -0.1812242994255389};
        Double[] c4RMDouble = new Double[] {c4RM.getEntry(0,0), c4RM.getEntry(0,1), c4RM.getEntry(1,1)};
        Assert.assertArrayEquals(c4Mat, c4RMDouble);

        // Matrix E
        RealMatrix e1RM = nodeMath.getEMatForNode(0);
        Double[] e1Mat = new Double[] {0.23952748745735447, -0.18127503919872656, -0.16434951438938683, 0.124380251933988};
        Double[] e1RMDouble = new Double[] {e1RM.getEntry(0,0), e1RM.getEntry(0,1), e1RM.getEntry(1,0), e1RM.getEntry(1,1)};
        Assert.assertArrayEquals(e1Mat, e1RMDouble);

        RealMatrix e2RM = nodeMath.getEMatForNode(1);
        Double[] e2Mat = new Double[] {0.23952748745735447, -0.18127503919872656, -0.16434951438938683, 0.124380251933988};
        Double[] e2RMDouble = new Double[] {e2RM.getEntry(0,0), e2RM.getEntry(0,1), e2RM.getEntry(1,0), e2RM.getEntry(1,1)};
        Assert.assertArrayEquals(e2Mat, e2RMDouble);

        RealMatrix e3RM = nodeMath.getEMatForNode(2);
        Double[] e3Mat = new Double[] {0.04666326455145753, -0.03531489040130742, -0.0320175632133628, 0.024230982265503798};
        Double[] e3RMDouble = new Double[] {e3RM.getEntry(0,0), e3RM.getEntry(0,1), e3RM.getEntry(1,0), e3RM.getEntry(1,1)};
        Assert.assertArrayEquals(e3Mat, e3RMDouble);

        RealMatrix e4RM = nodeMath.getEMatForNode(3);
        Double[] e4Mat = new Double[] {0.12398836438280014, -0.1666540537284611, -0.08507345812500677, 0.11434812323385568};
        Double[] e4RMDouble = new Double[] {e4RM.getEntry(0,0), e4RM.getEntry(0,1), e4RM.getEntry(1,0), e4RM.getEntry(1,1)};
        Assert.assertArrayEquals(e4Mat, e4RMDouble);

        // Vector b
        RealVector b1Vec = nodeMath.getbVecForNode(0);
        Double[] b1 = new Double[] {0.17091796790410735, 0.45694975625139694};
        Double[] b1RMDouble = new Double[] {b1Vec.getEntry(0), b1Vec.getEntry(1)};
        Assert.assertArrayEquals(b1, b1RMDouble);

        RealVector b2Vec = nodeMath.getbVecForNode(1);
        Double[] b2 = new Double[] {0.17091796790410735, 0.45694975625139694};
        Double[] b2RMDouble = new Double[] {b2Vec.getEntry(0), b2Vec.getEntry(1)};
        Assert.assertArrayEquals(b2, b2RMDouble);

        RealVector b3Vec = nodeMath.getbVecForNode(2);
        Double[] b3 = new Double[] {0.1922578250850458, 0.4407997494469062};
        Double[] b3RMDouble = new Double[] {b3Vec.getEntry(0), b3Vec.getEntry(1)};
        Assert.assertArrayEquals(b3, b3RMDouble);

        RealVector b4Vec = nodeMath.getbVecForNode(3);
        Double[] b4 = new Double[] {1.9305092441211686, 4.30118317706933};
        Double[] b4RMDouble = new Double[] {b4Vec.getEntry(0), b4Vec.getEntry(1)};
        Assert.assertArrayEquals(b4, b4RMDouble);

        // Vector d
        RealVector d1Vec = nodeMath.getdVecForNode(0);
        Double[] d1 = new Double[] {0.08879586049666498, -0.06092650619664548};
        Double[] d1RMDouble = new Double[] {d1Vec.getEntry(0), d1Vec.getEntry(1)};
        Assert.assertArrayEquals(d1, d1RMDouble);

        RealVector d2Vec = nodeMath.getdVecForNode(1);
        Double[] d2 = new Double[] {0.08879586049666498, -0.06092650619664548};
        Double[] d2RMDouble = new Double[] {d2Vec.getEntry(0), d2Vec.getEntry(1)};
        Assert.assertArrayEquals(d2, d2RMDouble);

        RealVector d3Vec = nodeMath.getdVecForNode(2);
        Double[] d3 = new Double[] {0.11038392303544745, -0.07573889798535485};
        Double[] d3RMDouble = new Double[] {d3Vec.getEntry(0), d3Vec.getEntry(1)};
        Assert.assertArrayEquals(d3, d3RMDouble);

        RealVector d4Vec = nodeMath.getdVecForNode(3);
        Double[] d4 = new Double[] {0.14214904428434422, -0.09753423969293895};
        Double[] d4RMDouble = new Double[] {d4Vec.getEntry(0), d4Vec.getEntry(1)};
        Assert.assertArrayEquals(d4, d4RMDouble);

        // double value f
        Assert.assertEquals(-3.24809910801962, nodeMath.getfForNode(0), EPSILON);
        Assert.assertEquals(-3.24809910801962, nodeMath.getfForNode(1), EPSILON);
        Assert.assertEquals( -4.87883483414817, nodeMath.getfForNode(2), EPSILON);
        Assert.assertEquals(-4.00043401352834, nodeMath.getfForNode(3), EPSILON);

        // Matrix L
        RealMatrix l1RM = nodeMath.getLMatForNode(0);
        Double[] l1Mat = new Double[] {-0.3757163513481878, 0.25779426672220895, -0.17688312929786404};
        Double[] l1RMDouble = new Double[] {l1RM.getEntry(0,0), l1RM.getEntry(0,1), l1RM.getEntry(1,1)};
        Assert.assertArrayEquals(l1Mat, l1RMDouble);

        RealMatrix l2RM = nodeMath.getLMatForNode(1);
        Double[] l2Mat = new Double[] {-0.3757163513481878, 0.25779426672220895, -0.17688312929786404};
        Double[] l2RMDouble = new Double[] {l2RM.getEntry(0,0), l2RM.getEntry(0,1), l2RM.getEntry(1,1)};
        Assert.assertArrayEquals(l2Mat, l2RMDouble);

        RealMatrix l3RM = nodeMath.getLMatForNode(2);
        Double[] l3Mat = new Double[] {-0.36977746373559073, 0.25371935362506814, -0.17408716516578368};
        Double[] l3RMDouble = new Double[] {l3RM.getEntry(0,0), l3RM.getEntry(0,1), l3RM.getEntry(1,1)};
        Assert.assertArrayEquals(l3Mat, l3RMDouble);

        //RealMatrix l4RM = nodeMath.getLMatForNode(3);
        //Double[] l4Mat = new Double[] {-0.746824771165506, 0.51242684261254, -0.351596892829451};
        //Double[] l4RMDouble = new Double[] {l4RM.getEntry(0,0), l4RM.getEntry(0,1), l4RM.getEntry(1,1)};
        //Assert.assertArrayEquals(l4Mat, l4RMDouble);

        RealMatrix l5RM = nodeMath.getLMatForNode(4);
        Double[] l5Mat = new Double[] {-0.7468247711655002, 0.5124268426125358, -0.3515968928294475};
        Double[] l5RMDouble = new Double[] {l5RM.getEntry(0,0), l5RM.getEntry(0,1), l5RM.getEntry(1,1)};
        Assert.assertArrayEquals(l5Mat, l5RMDouble);

        // Vector m
        RealVector m1Vec = nodeMath.getMVecForNode(0);
        Double[] m1 = new Double[] {0.37888042138556954, -0.2599648267685111};
        Double[] m1Double = new Double[] {m1Vec.getEntry(0), m1Vec.getEntry(1)};
        Assert.assertArrayEquals(m1, m1Double);

        RealVector m2Vec = nodeMath.getMVecForNode(1);
        Double[] m2 = new Double[] {0.3836914306467637, -0.26326584297153716};
        Double[] m2Double = new Double[] {m2Vec.getEntry(0), m2Vec.getEntry(1)};
        Assert.assertArrayEquals(m2, m2Double);

        RealVector m3Vec = nodeMath.getMVecForNode(2);
        Double[] m3 = new Double[] {0.16764381322398245, -0.11502723692442945};
        Double[] m3Double = new Double[] {m3Vec.getEntry(0), m3Vec.getEntry(1)};
        Assert.assertArrayEquals(m3, m3Double);

        //RealVector m4Vec = nodeMath.getMVecForNode(3);
        //Double[] m4 = new Double[] {0.170655395969879, -0.117093606302271};
        //Double[] m4Double = new Double[] {m4Vec.getEntry(0), m4Vec.getEntry(1)};
        //Assert.assertArrayEquals(m4, m4Double);

        RealVector m5Vec = nodeMath.getMVecForNode(4);
        Double[] m5 = new Double[] {0.3382992091938595, -0.23212084322670018};
        Double[] m5Double = new Double[] {m5Vec.getEntry(0), m5Vec.getEntry(1)};
        Assert.assertArrayEquals(m5, m5Double);

        Assert.assertEquals(-13.0101512466986, nodeMath.getRForNode(0), EPSILON);
        Assert.assertEquals(-13.6007534562081, nodeMath.getRForNode(1), EPSILON);
        Assert.assertEquals(-27.4541381503069, nodeMath.getRForNode(2), EPSILON);
        //Assert.assertEquals(-28.4013287142813, nodeMath.getRForNode(3), EPSILON);
        Assert.assertEquals(-55.8554668645882, nodeMath.getRForNode(4), EPSILON);


    }

}
