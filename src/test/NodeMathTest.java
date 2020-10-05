package test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.math.NodeMath;
import contraband.utils.PruneLikelihoodUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

/*
 * This class contains unit tests for NodeMath
 */
public class NodeMathTest {
    private int nTraits;
    private List<Double> data;
    private String spNames;
    private TreeParser tree;
    private String treeStr;
    private final KeyRealParameter traitValues = new KeyRealParameter();
    private final NodeMath nodeMath = new NodeMath();
    private double[] traitValuesArr;
    private RealParameter sigmasq;
    private RealParameter correlation;
    final static double EPSILON = 1e-7;
    /*
     * This test includes shrinkage estimations of trait correlations and operations of trait rate matrix
     */
    @Test
    public void testShrinkageOperations () {
        // tree
        treeStr = "((A:12.4420263,B:12.4420263):42.9258211,C:43.5702874);";
        spNames = "A B C";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 2;
        data = Arrays.asList(
                1.0, 2.0,
                3.0, 5.0,
                2.0, 4.0
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);
        traitValuesArr = new double[3 * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitValues, tree, nTraits, traitValuesArr);

        // populate a RealMatrix for trait value
        // the row number is the number of species
        // the column number is the number of traits
        RealMatrix traitRM = new Array2DRowRealMatrix(new double[3][nTraits]);
        PruneLikelihoodUtils.populateTraitValuesMatrix(traitValues, tree, nTraits, traitRM);

        // node math
        sigmasq = new RealParameter(new Double[]{0.1543038});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "shrinkage", true, "oneRateOnly", true);

        // (1) test the traditional unbiased estimation of trait correlations
        nodeMath.estimateCorrelations(traitRM);
        RealMatrix unBiasedEstimation = nodeMath.getUnbiasedRho();
        Assert.assertEquals(0.9819805060619656, unBiasedEstimation.getEntry(0,1), 0.0);

        // (2) test shrinkage estimation
        // the shrinkageEstimation is nTraits * nTraits real matrix with diagonal elements being 1
        // here we have two traits, there is 1 correlation
        nodeMath.populateShrinkageEstimation(0.25925925925926);
        RealMatrix shrinkageEstimation = nodeMath.getShrinkageRho();
        Assert.assertEquals(0.727392967453307, shrinkageEstimation.getEntry(0, 1), EPSILON);

        // (3) test the determinant of shrinkage matrix
        // the determinant is in real space
        double detShrinkageRho = nodeMath.getDetShrinkageRho();
        double detInvShrinkageRho = nodeMath.getDetInvShrinkageRho();
        Assert.assertEquals(0.470899470899473, detShrinkageRho, EPSILON);
        Assert.assertEquals(2.123595505617968, detInvShrinkageRho , EPSILON);

        // (4) test the trait rate matrix
        // When using shrinkage method, the trait rate matrix is used to transform the trait values
        nodeMath.populateTraitRateMatrix();
        RealMatrix traitRateRM = nodeMath.getTraitRateRealMatrix();
        Assert.assertEquals(0.1122394989713215, traitRateRM.getEntry(0, 1), EPSILON);

        // (5) test the determinant of trait rate matrix
        // When using shrinkage method, the determinant is calculated based on shrinkage rho matrix
        // the determinant is in log space
        nodeMath.populateInverseTraitRateMatrix();
        double detTraitRateMatrix = nodeMath.getTraitRateMatrixDeterminant();
        double detInvTraitRateMatrix = nodeMath.getTraitRateMatrixInverseDeterminant();
        Assert.assertEquals(-4.4907744304614, detTraitRateMatrix, EPSILON);
        Assert.assertEquals(4.4907744304614, detInvTraitRateMatrix, EPSILON);

        // (6) test the transformed trait values that are independent with other
        nodeMath.populateTransformedTraitValues(traitRM);
        double[] transformedTraitValues = nodeMath.getTransformedTraitValues();
        Assert.assertArrayEquals(new Double[] {2.5457261826187922, 4.721085948047915, 7.637178547856377, 10.45348269179349, 5.0914523652375845, 9.44217189609583}, ArrayUtils.toObject(transformedTraitValues));
    }

    /*
     * This test takes input trait correlations and populate trait rate matrix
     * Tests include operations of trait rate matrix
     *
     * Note: all traits share the same evolutionary rate
     */
    @Test
    public void testNonShrinkageOperationsOneRateOnly () {
        // tree
        treeStr = "((A:23.0058179,B:23.0058179):14.350951,C:37.3567689);";
        spNames = "A B C";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 2;
        data =  Arrays.asList(
                -2.62762948691895, -1.56292164859448,
                -1.50846427625826, -1.59482814741543,
                -0.226074849617958, -2.11000367246907
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);
        traitValuesArr = new double[3 * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitValues, tree, nTraits, traitValuesArr);

        // BM model parameters
        sigmasq = new RealParameter(new Double[]{0.3145740});
        correlation = new RealParameter(new Double[]{-0.632620487603683});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation, "oneRateOnly", true);

        // (1) test the trait rate matrix
        nodeMath.populateTraitRateMatrixForOneRate();
        double [] traitRateMatrix = nodeMath.getTraitRateMatrix();
        Assert.assertArrayEquals(new Double[] {0.314574, -0.19900595726744097, -0.19900595726744097, 0.314574}, ArrayUtils.toObject(traitRateMatrix));

        // (2) test the determinant of trait rate matrix
        nodeMath.performMatrixOperations();
        Assert.assertEquals(-2.824245359292377, nodeMath.getTraitRateMatrixDeterminant(), EPSILON); // the determinant is in log space
        Assert.assertEquals(16.84822583043347, nodeMath.getTraitRateMatrixInverseDeterminant(), EPSILON); // the determinant is in real space
    }

    /*
     * This test takes input trait correlations and populate trait rate matrix
     * Tests include operations of trait rate matrix
     *
     * Note: every trait has its own evolutionary rate
     */
    @Test
    public void testNonShrinkageOperationsMultipleRates () {
        // tree
        treeStr = "((A:23.0058179,B:23.0058179):14.350951,C:37.3567689);";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 2;
        data = Arrays.asList(
                -2.62762948691895, -0.764018322006132,
                -1.50846427625826, -1.02686498716963,
                -0.226074849617958, -1.73165056392106
        );
        spNames = "A B C";
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);
        traitValuesArr = new double[3 * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitValues, tree, nTraits, traitValuesArr);

        // BM model parameters
        sigmasq = new RealParameter(new Double[] {0.3, 0.2});
        correlation = new RealParameter(new Double[] {-0.720107524122507});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation);

        // (1) test trait rate matrix
        nodeMath.populateTraitRateMatrixForMultipleRates();
        double [] traitRateMatrix = nodeMath.getTraitRateMatrix();
        Assert.assertArrayEquals(new Double[] { 0.3, -0.17638959940390708, -0.17638959940390708, 0.2 }, ArrayUtils.toObject(traitRateMatrix));

        // (2) test the determinant of trait rate matrix
        nodeMath.performMatrixOperations();
        Assert.assertEquals(-3.544373678152544, nodeMath.getTraitRateMatrixDeterminant(), 0.0); // the determinant is in log space
        Assert.assertEquals(34.61799654333531, nodeMath.getTraitRateMatrixInverseDeterminant(), 0.0); // the determinant is in real space
    }
}
