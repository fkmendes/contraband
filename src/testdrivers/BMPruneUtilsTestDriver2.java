package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.prunelikelihood.BMPruneLikelihood;
import contraband.utils.PruneLikelihoodUtils;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

public class BMPruneUtilsTestDriver2 {
    public static void main(String[] args) {
        int nTraits = 2;
        int nSpecies = 3;
        List<Double> data = Arrays.asList(
                -2.62762948691895, -0.764018322006132,
                -1.50846427625826, -1.02686498716963,
                -0.226074849617958, -1.73165056392106
        );
        String spNames = "A B C";

        KeyRealParameter traitsValues = new KeyRealParameter();
        traitsValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        String treeStr = "((A:23.0058179,B:23.0058179):14.350951,C:37.3567689);";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        double[] traitValuesArrayList = new double[nSpecies * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitsValues, tree, nTraits, traitValuesArrayList);

        RateCategoryClockModel  pcmc = new RateCategoryClockModel();
        IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
        RealParameter colorValues = new RealParameter(new Double[] {1.0});
        pcmc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        NodeMath nodeMath = new NodeMath();
        RealParameter sigmasq = new RealParameter(new Double[] {0.3, 0.2});
        RealParameter correlation = new RealParameter(new Double[] {-0.720107524122507});
        RealParameter rootValues = new RealParameter(new Double[] {-1.31465955080609, -1.2374605274288});
        nodeMath.initByName("traits", traitsValues, "sigmasq", sigmasq, "correlation", correlation, "rootValues", rootValues);

        nodeMath.performMatrixOperations();
        System.out.println(Arrays.toString(nodeMath.getTraitRateMatrix()));

        System.out.println("Print detTraitRateMatrix = " + nodeMath.getTraitRateMatrixDeterminant());
        // expected: 0.0288867092221292

        System.out.println("Print detTraitRateMatrix = " + nodeMath.getTraitRateMatrixInverseDeterminant());
        // expected: 34.6179965433353

        BMPruneLikelihood pcm = new BMPruneLikelihood();
        pcm.pruneNode(tree.getRoot(), nTraits, traitValuesArrayList, pcmc, nodeMath, false);

        int rootIdx = tree.getRoot().getNr();

        double l0 = nodeMath.getLForNode(rootIdx);
        System.out.println("Print L at root = " + l0);
        // expected:

        // get the mVec at the root
        double[] m0 = nodeMath.getMVecForNode(rootIdx);
        System.out.println("Print mVec at root = " + Arrays.toString(m0));
        // expected:

        double r0 = nodeMath.getRForNode(rootIdx);
        System.out.println("Print r at root =" + r0);
        // expected:

        // lMat.E * t(mE) %*% Rho.inverse %*% mE + t(mE) %*% mVec.E + r.E
        double logP = l0 * MatrixUtilsContra.tVecDotMatrixDotVec(nodeMath.getRootValuesArr(), nodeMath.getTraitRateMatrixInverse(), nTraits) + MatrixUtilsContra.vectorDotMultiply(nodeMath.getRootValuesArr(), m0) + r0;

        System.out.println("Print log likelihood = " + logP);
        // expected: -11.0563970153267
    }
}
