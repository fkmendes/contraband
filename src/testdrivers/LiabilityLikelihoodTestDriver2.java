package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.NodeMath;
import contraband.prunelikelihood.LiabilityLikelihood;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class LiabilityLikelihoodTestDriver2 {
    public static void main(String[] args) {

        // tree
        String treeStr = "((t1:0.1252127116,t2:0.1252127116):0.1805853832,(t3:0.03846640117,t4:0.03846640117):0.2673316937);";
        String spNames = "t1 t2 t3 t4";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        int nTraits = 3;
        List<Double> data = Arrays.asList(
                -0.844322389522931, -0.19057708979688, 1.14244560543034,
                -1.82272722113192, -0.451751783026192, 1.4540953632918,
                -0.402021030584138, -0.399107913371344, -0.546709817857608,
                -0.653736939856334, -0.567529184656187, -0.862192546948912
        );
        KeyRealParameter traitValues = new KeyRealParameter();
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        NodeMath nodeMath = new NodeMath();
        // BM model parameters
        RealParameter inverseMatrix = new RealParameter(new Double[] {
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0
        });
        // all traits share one rate
        RealParameter sigmasq = new RealParameter(new Double[] {1.0});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "inverseMatrix", inverseMatrix, "oneRateOnly", true);

        // branch rate model
        IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
        RealParameter colorValues = new RealParameter(new Double[] {1.0});
        RateCategoryClockModel lsc = new RateCategoryClockModel();
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // non-transformed likelihood
        LiabilityLikelihood lik1 = new LiabilityLikelihood();
        lik1.initByName("traits", traitValues,
                "tree", tree, "nodeMath", nodeMath, "branchRateModel", lsc);

        double logP1 = lik1.calculateLogP();
        System.out.println("Log likelihood 1 = "+ logP1); //-9.01503846914647

        // transformed likelihood
        NodeMath nodeMath2 = new NodeMath();
        nodeMath2.initByName("traits", traitValues, "sigmasq", sigmasq, "inverseMatrix", inverseMatrix, "oneRateOnly", true, "transform", true);
        LiabilityLikelihood lik2 = new LiabilityLikelihood();
        lik2.initByName("traits", traitValues,
                "tree", tree, "nodeMath", nodeMath2, "branchRateModel", lsc);
        double logP2 = lik2.calculateLogP();
        System.out.println("Log likelihood 2 = "+ logP2); //-9.01503846914647
    }
}
