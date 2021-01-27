package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.GeneralNodeMath;
import contraband.prunelikelihood.BMLikelihood;
import contraband.prunelikelihood.MorphologicalData;
import contraband.prunelikelihood.SigmaMatrix;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class BMLikelihoodTestDriver {

    public static void main(String[] args) {
        // tree
        String treeStr = "(((t3:1.209461463,t4:1.209461463):0.6659547705,(t9:0.841016425,t10:0.841016425):1.034399809):1.561956365,(((t2:1.602817551,(t5:1.164343725,t6:1.164343725):0.4384738261):0.6643605462,t1:2.267178098):0.7120187616,(t7:1.115655119,t8:1.115655119):1.86354174):0.458175739);";
        String spNames = "t3 t4 t9 t10 t2 t5 t6 t1 t7 t8";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        int nTraits = 2;
        KeyRealParameter traitValues = new KeyRealParameter();
        List<Double> data =  Arrays.asList(
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

        MorphologicalData morphData = new MorphologicalData();
        morphData.initByName("traits", traitValues, "tree", tree);

        // branch rate model
        RateCategoryClockModel lsc = new RateCategoryClockModel();
        IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
        RealParameter colorValues = new RealParameter(new Double[] {1.0});
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // BM model parameters
        RealParameter sigmasq = new RealParameter(new Double[]{0.93818618860063, 1.74928729364035});
        RealParameter correlation = new RealParameter(new Double[]{-0.117852694436286});
        SigmaMatrix sigmaMatrix = new SigmaMatrix();
        sigmaMatrix.initByName("sigmasq", sigmasq, "correlation", correlation, "trait", morphData);

        GeneralNodeMath nodeMath = new GeneralNodeMath();
        nodeMath.initByName("trait", morphData, "rateMatrix", sigmaMatrix, "tree", tree);

        // prune likelihood
        BMLikelihood pcm = new BMLikelihood();
        pcm.initByName("nodeMath", nodeMath, "tree", tree, "trait", morphData, "branchRateModel", lsc);
        double logP = pcm.calculateLogP();
        System.out.println("Likelihood = " + logP); // -39.588461762377


    }
}
