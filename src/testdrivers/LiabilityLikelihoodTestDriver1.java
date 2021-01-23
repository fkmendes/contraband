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

public class LiabilityLikelihoodTestDriver1 {
    public static void main(String[] args) {
        // tree
        String treeStr = "((t4:0.294791968,(t7:0.2380450868,t8:0.2380450868):0.05674688124):1.269758843,(((t5:0.2457411247,t6:0.2457411247):0.4091478892,(t2:0.3449767887,t3:0.3449767887):0.3099122252):0.771197872,(t1:1.012806474,(t9:0.06711245421,t10:0.06711245421):0.9456940195):0.4132804122):0.1384639249);";
        String spNames = "t4 t7 t8 t5 t6 t2 t3 t1 t9 t10";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        int nTraits = 4;
        List<Double> data = Arrays.asList(
                -6.27499208443443, -1.22468097213209, 5.81907241368362, -3.92206107638759,
                -4.94850658278194, -0.453876476004197, 4.30633969774666, -3.3563174816445,
                -6.02002559472886, -2.09468363352179, 4.46574877202188, -4.26072495377516,
                -5.0798244799288, -1.72613398635412, 5.76821474836136, 0.707038879576404,
                -5.50430467741627, -3.14168739896938, 4.56177348682651, 3.9073600922006,
                -4.12305027861692, -3.23866814100161, 3.9661783662543, 2.52600581804315,
                -5.3875721133069, -5.37340326805192, 5.70660547120509, 0.111472504390871,
                -3.21631588401429, -4.52403176018136, 7.3156952337464, -1.77999035455415,
                -5.02598308954637, -4.26615711744045, 1.1431561920545, -0.861746816351694,
                -6.68503249988375, -5.39752125164611, 1.27697631062818, -0.853488211694949
        );
        KeyRealParameter traitValues = new KeyRealParameter();
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        NodeMath nodeMath = new NodeMath();
        // BM model parameters
        RealParameter inverseMatrix = new RealParameter(new Double[] {
                1.15967096041944, -0.22059864425817, 0.0349183919101759, 0.457034581560665,
                -0.22059864425817, 1.19212258045853, 0.210910014702716, -0.359004307577557,
                0.0349183919101759, 0.210910014702716, 1.27408118183747, 0.496761070565367,
                0.457034581560665, -0.359004307577557, 0.496761070565367, 1.47635603144872
        });
        // all traits share one rate
        RealParameter sigmasq = new RealParameter(new Double[] {4.91917941793269});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "inverseMatrix", inverseMatrix, "oneRateOnly", true, "transform", true);

        // branch rate model
        IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
        RealParameter colorValues = new RealParameter(new Double[] {1.0});
        RateCategoryClockModel lsc = new RateCategoryClockModel();
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // transformed likelihood
        LiabilityLikelihood lik = new LiabilityLikelihood();
        lik.initByName("traits", traitValues,
                "tree", tree, "nodeMath", nodeMath, "branchRateModel", lsc);
        double logP = lik.calculateLogP();
        System.out.println("Log likelihood1 = "+ logP); // -75.80052159761348

        // non-transformed likelihood
        NodeMath nodeMath2 = new NodeMath();
        nodeMath2.initByName("traits", traitValues, "sigmasq", sigmasq, "inverseMatrix", inverseMatrix, "oneRateOnly", true);
        // transformed likelihood
        LiabilityLikelihood lik2 = new LiabilityLikelihood();
        lik2.initByName("traits", traitValues,
                "tree", tree, "nodeMath", nodeMath2, "branchRateModel", lsc);
        double logP2 = lik2.calculateLogP();
        System.out.println("Log likelihood2 = "+ logP2); // -75.80052159761351
    }
}
