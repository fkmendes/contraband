package contraband.test;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.NodeMath;
import contraband.prunelikelihood.BMPruneLikelihood;
import contraband.prunelikelihood.BMPruneShrinkageLikelihood;
import contraband.utils.PruneLikelihoodUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;
import java.util.Arrays;
import java.util.List;

/*
 * This class contains unit tests for PruneLikelihoodProcess,
 * which prunes the tree and calculates the likelihood for continuous traits.
 */

public class PruneLikelihoodProcessTest {

    private TreeParser tree;
    private String treeStr;
    private Integer nTraits;
    private String spNames;
    private List<Double> data;
    private final RealParameter traitValues = new RealParameter();

    /*
     * Tree with 3 species and each species has 2 traits
     * Traits share one evolutionary rate
     * Root values are estimated
     * Without using shrinkage method, i.e. trait correlation is taken from input
     */
    @Test
    public void testOneRateOnlyWithoutShrinkage () {
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

        // populate an array for trait values, which includes elements of a 3 * 2 matrix
        double[] traitValuesArray = new double[3 * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitValues, tree, nTraits, traitValuesArray);

        // BM model parameters
        RealParameter sigmasq = new RealParameter(new Double[]{0.3145740}); // shared evolutionary rate
        RealParameter correlation = new RealParameter(new Double[]{-0.632620487603683}); // trait correlation
        NodeMath nodeMath = new NodeMath();
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation,"oneRateOnly", true);

        // the operations include calculating the determinant of trait rate matrix and the determinant of inverse trait rate matrix
        // which is necessary for likelihood calculation
        nodeMath.performMatrixOperations();

        // branch rate model
        RateCategoryClockModel pcmc = new RateCategoryClockModel();
        IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
        RealParameter colorValues = new RealParameter(new Double[] {1.0});
        pcmc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // run pruning algorithm
        BMPruneLikelihood pcm = new BMPruneLikelihood();
        pcm.pruneNode(tree.getRoot(), nTraits, traitValuesArray, pcmc, nodeMath, false);

        int rootIdx = tree.getRoot().getNr();
        // l0: l value at the root
        Assert.assertEquals(-0.032723927183444676, nodeMath.getLForNode(rootIdx), 0.0);
        // m0: m vector at the root
        Assert.assertArrayEquals(new Double[] {-0.8501607370417381, -0.9115145062074426}, ArrayUtils.toObject(nodeMath.getMVecForNode(rootIdx)));
        // r0: r value at the root
        Assert.assertEquals(-13.528327328718111, nodeMath.getRForNode(rootIdx), 0.0);

        // if root value are not specified as an input of NodeMath
        // the trait values at the root are maximum likelihood estimations
        nodeMath.populateRootValuesVec(rootIdx);
        double [] rootValuesArray = nodeMath.getRootValuesArr();
        Assert.assertArrayEquals(new Double[] {-1.3146595508060854, -1.7961125556621191}, ArrayUtils.toObject(rootValuesArray));
    }

    /*
     * Tree with 3 species and each species has 2 traits
     * Traits have own evolutionary rates
     * Root values are specified
     * Without using shrinkage method, i.e. trait correlation is taken from input
     */
    @Test
    public void testMultipleRatesWithoutShrinkage () {
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

        // populate an array for trait values, which includes elements of a 3 * 2 matrix
        double[] traitValuesArray = new double[3 * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitValues, tree, nTraits, traitValuesArray);

        // branch rate model
        RateCategoryClockModel  pcmc = new RateCategoryClockModel();
        IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
        RealParameter colorValues = new RealParameter(new Double[] {1.0});
        pcmc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        NodeMath nodeMath = new NodeMath();
        RealParameter sigmasq = new RealParameter(new Double[] {0.3, 0.2}); // evolutionary rates
        RealParameter correlation = new RealParameter(new Double[] {-0.720107524122507}); // trait correlation between the two traits
        RealParameter rootValues = new RealParameter(new Double[] {-1.31465955080609, -1.2374605274288}); // root values are specified as an input
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation, "rootValues", rootValues,"oneRateOnly", false);

        // the operations include calculating the determinant of trait rate matrix and the determinant of inverse trait rate matrix
        // which is necessary for likelihood calculation
        nodeMath.performMatrixOperations();

        // run pruning algorithm
        BMPruneLikelihood pcm = new BMPruneLikelihood();
        pcm.pruneNode(tree.getRoot(), nTraits, traitValuesArray, pcmc, nodeMath, false);

        int rootIdx = tree.getRoot().getNr();
        // l0: l value at the root
        Assert.assertEquals(-0.032723927183444676, nodeMath.getLForNode(rootIdx), 0.0);
        // m0: m vector at the root
        Assert.assertArrayEquals(new Double[] { -1.0902581683947343, -1.3664966897695918 }, ArrayUtils.toObject(nodeMath.getMVecForNode(rootIdx)));
        // r0: r value at the root
        Assert.assertEquals(-12.618549029514906, nodeMath.getRForNode(rootIdx), 0.0);
    }

    /*
     * Tree with 3 species and each species has 2 traits
     * Traits share one evolutionary rate
     * Root values are estimated
     * Shrinkage method is used, i.e. trait correlation is estimated
     */
    @Test
    public void testOneRateOnlyWithShrinkage () {
        // tree
        treeStr = "((A:12.4420263,B:12.4420263):42.9258211,C:43.5702874);";
        spNames = "A B C";
        tree = new TreeParser(treeStr, false, false, true, 0);

        //trait values
        nTraits = 2;
        data = Arrays.asList(
                1.0, 2.0,
                3.0, 5.0,
                2.0, 4.0
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // populate a 3 * 2 real matrix for trait values
        // which will be used in unbiased estimation of trait correlation
        RealMatrix traitRM = new Array2DRowRealMatrix(new double[3][nTraits]);
        PruneLikelihoodUtils.populateTraitValuesMatrix(traitValues, tree, nTraits, traitRM);

        // node math
        RealParameter sigmasq = new RealParameter(new Double[]{0.1543038}); // shared evolutionary rate
        NodeMath nodeMath = new NodeMath();
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "shrinkage", true, "oneRateOnly", true);

        // unbiased estimation of trait correlation
        nodeMath.estimateCorrelations(traitRM);
        // shrinkage estimation given the shrinkage parameter
        nodeMath.populateShrinkageEstimation(0.25925925925926);

        // use the estimated correlation and the specified evolutionary rate to populate the trait rate matrix
        // then calculate the determinant of trait rate matrix and the determinant of inverse trait rate matrix
        nodeMath.populateTraitRateMatrix();
        nodeMath.populateInverseTraitRateMatrix();

        // use the trait rate matrix to transform the trait values
        // so that the traits used in the likelihood calculation are independent from each other
        nodeMath.populateTransformedTraitValues(traitRM);

        // branch rate model
        RateCategoryClockModel pcmc = new RateCategoryClockModel();
        IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
        RealParameter colorValues = new RealParameter(new Double[] {1.0});
        pcmc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // PruneLikelihood with shrinkage
        BMPruneShrinkageLikelihood pcm = new BMPruneShrinkageLikelihood();
        // ignore the population variance
        pcm.setPopSE(false);
        // run pruning algorithm
        pcm.pruneNode(tree.getRoot(), nTraits, nodeMath.getTransformedTraitValues(), pcmc, nodeMath, false);

        int rootIdx = tree.getRoot().getNr();
        // l0: l value at the root
        Assert.assertEquals(-0.02164930565494149, nodeMath.getLForNode(rootIdx), 0.0);
        // m0: m vector at the root
        Assert.assertArrayEquals(new Double[] {0.22045281696520652, 0.37109117994225865}, ArrayUtils.toObject(nodeMath.getMVecForNode(rootIdx)));
        // r0: r value at the root
        Assert.assertEquals(-13.012014942429309, nodeMath.getRForNode(rootIdx), 0.0);

        // if root value are not specified as an input of NodeMath
        // the trait values at the root are maximum likelihood estimations
        nodeMath.populateRootValuesVec(rootIdx);
        double [] rootValuesArray = nodeMath.getRootValuesArr();
        Assert.assertArrayEquals(new Double[] {5.0914523652375845, 8.570509970548557}, ArrayUtils.toObject(rootValuesArray));
    }
}
