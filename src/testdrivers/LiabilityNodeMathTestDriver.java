package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.LiabilityNodeMath;
import contraband.prunelikelihood.BMPruneLikelihood;
import contraband.prunelikelihood.LiabilityLikelihood;
import org.junit.Assert;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

public class LiabilityNodeMathTestDriver {

    public static void main(String[] args) {
        // tree
        String treeStr = "((t4:0.294791968,(t7:0.2380450868,t8:0.2380450868):0.05674688124):1.269758843,(((t5:0.2457411247,t6:0.2457411247):0.4091478892,(t2:0.3449767887,t3:0.3449767887):0.3099122252):0.771197872,(t1:1.012806474,(t9:0.06711245421,t10:0.06711245421):0.9456940195):0.4132804122):0.1384639249);";
        String spNames = "t4 t7 t8 t5 t6 t2 t3 t1 t9 t10";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        int nTraits = 4;
        List<Double> data = Arrays.asList(
                -6.42193809335608, -1.65989756061141, 5.15759789266981, -6.56062921109306,
                -5.03724385041298, 0.93007037382525, 4.41105696963357, -5.84790423473936,
                -6.15578318250413, -3.89561377661385, 3.84425761908253, -4.58575554859401,
                -5.17432423006894, -2.53719420891046, 5.03090939175046, -0.79379415322876,
                -5.61743143702448, -6.50421680550017, 3.57379452807292, 6.3680803136263,
                -4.1755649321327, -6.22257127540476, 3.20155065335743, 4.7549882614031,
                -5.49557642294037, -12.4500805740433, 3.6241593165236, 3.17971650657134,
                -3.2290412841014, -9.31768106789294, 5.31489864149079, -2.96128941493352,
                -5.11812017566976, -9.33547489795037, 0.63612585434188, 5.31629275763256,
                -6.84997187279757, -13.0239682050792, 0.189981169418762, 7.38483491142637
        );
        KeyRealParameter traitValues = new KeyRealParameter();
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        LiabilityNodeMath nodeMath = new LiabilityNodeMath();
        // BM model parameters
        RealParameter inverseMatrix = new RealParameter(new Double[] {
                0.431864494796763, 0.269799734818774, 0.265743456802657
        });
        // all traits share one rate
        RealParameter sigmasq = new RealParameter(new Double[] {1.0});

        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "inverseMatrix", inverseMatrix);

        // branch rate model

        IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
        RealParameter colorValues = new RealParameter(new Double[] {1.0});
        RateCategoryClockModel lsc = new RateCategoryClockModel();
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // prune likelihood
        LiabilityLikelihood liabilityLikelihood = new LiabilityLikelihood();
        liabilityLikelihood.initByName("traits", traitValues,
                "tree", tree, "nodeMath", nodeMath, "branchRateModel", lsc);

        double lik = liabilityLikelihood.calculateLogP();
        System.out.println("Log likelihood"+ lik);
        
    }
}
