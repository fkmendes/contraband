package testdrivers;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeParser;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.clock.RateCategoryClockModel;
import contraband.prunelikelihood.BMPruneLikelihood;
import contraband.utils.PruneLikelihoodUtils;
import java.util.Arrays;
import java.util.List;

public class BMPruneUtilsTestDriver1 {
    public static void main(String[] args) {
        int nTraits = 2;
        int nSpecies = 3;
        List<Double> data = Arrays.asList(
                -2.62762948691895, -1.56292164859448,
                -1.50846427625826, -1.59482814741543,
                -0.226074849617958, -2.11000367246907
        );
        String spNames = "A B C";

        RealParameter traitsValues = new RealParameter();
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
        RealParameter sigmasq = new RealParameter(new Double[] {0.3145740});
        RealParameter correlation = new RealParameter(new Double[] {-0.632620487603683});
        RealParameter rootValues = new RealParameter(new Double[] {-1.31465955080609, -1.79611255566212});
        nodeMath.initByName("traits", traitsValues, "sigmasq", sigmasq, "correlation", correlation, "rootValues", rootValues, "oneRateOnly", true);

        nodeMath.performMatrixOperations();

        System.out.println("Print detTraitRateMatrix = " + nodeMath.getTraitRateMatrixDeterminant());
        // expected: 0.599791318664079

        System.out.println("Print detTraitRateMatrix = " + nodeMath.getTraitRateMatrixInverseDeterminant());
        // expected: 1.66724653872502

        BMPruneLikelihood pcm = new BMPruneLikelihood();
        pcm.pruneNode(tree.getRoot(), nTraits, traitValuesArrayList, pcmc, nodeMath, false);

        int rootIdx = tree.getRoot().getNr();

        double l0 = nodeMath.getLForNode(rootIdx);
        System.out.println("Print L at root = " + l0);
        // expected: -0.032723927183444676

        // get the mVec at the root
        double[] m0 = nodeMath.getMVecForNode(rootIdx);
        System.out.println("Print mVec at root = " + Arrays.toString(m0));
        // expected: -0.850160737041739, -0.911514506207443

        double r0 = nodeMath.getRForNode(rootIdx);
        System.out.println("Print r at root =" + r0);
        // expected: -13.5283273287181

        // lMat.E * t(mE) %*% Rho.inverse %*% mE + t(mE) %*% mVec.E + r.E
        double logP = l0 * MatrixUtilsContra.tVecDotMatrixDotVec(nodeMath.getRootValuesArr(), nodeMath.getTraitRateMatrixInverse(), nTraits) + MatrixUtilsContra.vectorDotMultiply(nodeMath.getRootValuesArr(), m0) + r0;

        System.out.println("Print log likelihood = " + logP);
        // expected: -12.1509000377483
    }
}
