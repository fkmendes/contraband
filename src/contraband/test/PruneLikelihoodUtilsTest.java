package contraband.test;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeParser;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.utils.PruneLikelihoodUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;
import java.util.Arrays;
import java.util.List;

/*
 * This class contains unit tests for PruneLikelihoodUtils,
 * which populates continuous trait values and parameters for pruning likelihood calculation.
 */

public class PruneLikelihoodUtilsTest {

    private TreeParser tree;
    private String treeStr;
    private Integer nTraits;
    private String spNames;
    private List<Double> data;
    private final RealParameter traitValues = new RealParameter();
    final static double EPSILON = 1e-7;
    /*
     * (1) Populates trait values in an array and a real matrix
     */
    @Test
    public void testPopulateTraitValues (){
        // tree
        treeStr = "((sp3:0.1201847336,sp4:0.1201847336):0.8798152664,(sp1:0.9098543416,sp2:0.9098543416):0.09014565845):0.0;";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // four species
        // each species has two traits
        nTraits = 2;
        spNames = "sp3 sp4 sp1 sp2";
        int nSpecies = 4;
        // the data is given in accordance with the species names
        data =  Arrays.asList(
                0.983714690867666, -7.54729477473779,
                -7.86424514338822, -2.97908131550921,
                7.23079460908758, -0.780647498381348,
                -1.39605330265115, 3.72028693114977
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // (1) populates the trait values in an array
        double[] traitValuesArr = new double[nSpecies * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitValues,tree, nTraits, traitValuesArr);
        // the expected the array gives the trait values in accordance with the node number
        Assert.assertArrayEquals(new Double[] {7.23079460908758, -0.780647498381348, -1.39605330265115, 3.72028693114977, 0.983714690867666, -7.54729477473779, -7.86424514338822, -2.97908131550921}, ArrayUtils.toObject(traitValuesArr));

        // (2) populates the trait values in a real matrix
        // the column number is equal to traits number, i.e. 2
        // the row number is equal to species number, i.e. 4
        RealMatrix traitValuesRM = new Array2DRowRealMatrix(new double[nSpecies][nTraits]);
        PruneLikelihoodUtils.populateTraitValuesMatrix(traitValues, tree, nTraits, traitValuesRM);
        // sp1
        Assert.assertArrayEquals(new Double[] {7.23079460908758, -0.780647498381348}, ArrayUtils.toObject(traitValuesRM.getRow(0)));
        // sp2
        Assert.assertArrayEquals(new Double[] {-1.39605330265115, 3.72028693114977}, ArrayUtils.toObject(traitValuesRM.getRow(1)));
        // sp3
        Assert.assertArrayEquals(new Double[] {0.983714690867666, -7.54729477473779}, ArrayUtils.toObject(traitValuesRM.getRow(2)));
        // sp4
        Assert.assertArrayEquals(new Double[] {-7.86424514338822, -2.97908131550921}, ArrayUtils.toObject(traitValuesRM.getRow(3)));
    }

    /*
     * (2) Prune a tree with 3 species and 2 traits without shrinkage method
     */
    @Test
    public void testPopulatePCMParams () {
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
        double[] traitValuesArray = new double[3 * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitValues, tree, nTraits, traitValuesArray);

        // BM model parameters
        // When shrinkage method is not used, correlations among traits are taken from input
        RealParameter sigmasq = new RealParameter(new Double[]{0.3145740});
        RealParameter correlation = new RealParameter(new Double[]{-0.632620487603683});
        RealParameter rootValues = new RealParameter(new Double[] {-1.31465955080609, -1.79611255566212});
        NodeMath nodeMath = new NodeMath();
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation, "rootValues", rootValues,"oneRateOnly", true);

        // populate trait rate matrix and calculate determinant
        nodeMath.performMatrixOperations();

        // initialize the mVector for the parent of A and B
        double[] mVecD = new double [2];

        // the first node, i.e. tip A
        int tipIdx = 0;
        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(tipIdx).getLength(), nTraits, tipIdx);
        Assert.assertEquals(-0.021733632865102357, nodeMath.getAForNode(tipIdx), 0.0);
        Assert.assertEquals(-0.021733632865102357, nodeMath.getCForNode(tipIdx), 0.0);
        Assert.assertEquals(0.043467265730204714, nodeMath.getEForNode(tipIdx), 0.0);
        Assert.assertEquals(-3.561501522879213, nodeMath.getfForNode(tipIdx), EPSILON);

        PruneLikelihoodUtils.populateLmrForTip(nodeMath, traitValuesArray, nTraits, tipIdx);
        Assert.assertEquals(-0.021733632865102357,nodeMath.getLForNode(tipIdx), 0.0);
        Assert.assertArrayEquals(new Double[] {-0.8331278807826273, -0.7430154496439219}, ArrayUtils.toObject(nodeMath.getMVecForNode(tipIdx)));
        Assert.assertEquals(-5.2367146815829,nodeMath.getRForNode(tipIdx), EPSILON);

        MatrixUtilsContra.vectorAdd(mVecD, nodeMath.getMVecForNode(tipIdx), mVecD);

        // the second node, i.e. tip B
        tipIdx = 1;
        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(tipIdx).getLength(), nTraits, tipIdx);
        Assert.assertEquals(-0.02173363286510235, nodeMath.getAForNode(tipIdx), EPSILON);
        Assert.assertEquals(-0.02173363286510235, nodeMath.getCForNode(tipIdx), EPSILON);
        Assert.assertEquals(0.04346726573020471, nodeMath.getEForNode(tipIdx), EPSILON);
        Assert.assertEquals(-3.561501522879213, nodeMath.getfForNode(tipIdx), EPSILON);

        PruneLikelihoodUtils.populateLmrForTip(nodeMath, traitValuesArray, nTraits, tipIdx);
        Assert.assertEquals(-0.021733632865102357,nodeMath.getLForNode(tipIdx), EPSILON);
        Assert.assertArrayEquals(new Double[] {-0.5799479302276102, -0.5872574081072599}, ArrayUtils.toObject(nodeMath.getMVecForNode(tipIdx)));
        Assert.assertEquals(-4.467204212412191,nodeMath.getRForNode(tipIdx), EPSILON);

        // the mVector is added up to their parent
        MatrixUtilsContra.vectorAdd(mVecD, nodeMath.getMVecForNode(tipIdx), mVecD);

        // the third node, i.e. internal node that is the parent of A and B
        int internalIdx = 3;
        // first we need to set the L m r coming from its child nodes
        nodeMath.setLForNode(internalIdx , nodeMath.getLForNode(0) + nodeMath.getLForNode(1));
        nodeMath.setMVecForNode(internalIdx, mVecD);
        nodeMath.setRForNode(internalIdx , nodeMath.getRForNode(0) + nodeMath.getRForNode(1));
        Assert.assertEquals( -0.04346726573020471, nodeMath.getLForNode(internalIdx), EPSILON);
        Assert.assertArrayEquals(new Double[] {-1.4130758110102375, -1.3302728577511818}, ArrayUtils.toObject(nodeMath.getMVecForNode(internalIdx)));
        Assert.assertEquals(-9.703918893995091, nodeMath.getRForNode(internalIdx), 0.0);

        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(internalIdx).getLength(), nTraits, internalIdx);
        Assert.assertEquals(-0.03484089660678236, nodeMath.getAForNode(internalIdx), EPSILON);
        Assert.assertEquals(-0.03484089660678236, nodeMath.getCForNode(internalIdx), EPSILON);
        Assert.assertEquals(0.06968179321356473, nodeMath.getEForNode(internalIdx), EPSILON);
        Assert.assertEquals(-3.089570598549913, nodeMath.getfForNode(internalIdx), EPSILON);

        PruneLikelihoodUtils.populateLmrForInternalNode(nodeMath, nTraits,internalIdx);
        Assert.assertEquals(-0.0193394719769881, nodeMath.getLForNode(internalIdx), EPSILON);
        Assert.assertArrayEquals(new Double[] {-0.6287062134990086, -0.5918654928494784}, ArrayUtils.toObject(nodeMath.getMVecForNode(internalIdx)));
        Assert.assertEquals(-9.119795872783662, nodeMath.getRForNode(internalIdx), 0.0);
    }

    /*
     * (3) Prune a tree with 3 species and 2 traits with shrinkage method
     */
    @Test
    public void testPopulatePCMParamsWithShrinkage () {
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

        RealMatrix traitValuesRM = new Array2DRowRealMatrix(new double[3][nTraits]);
        PruneLikelihoodUtils.populateTraitValuesMatrix(traitValues, tree, nTraits, traitValuesRM);

        // BM model parameters
        // When shrinkage method is used, correlations among traits are estimated in NodeMath
        RealParameter sigmasq = new RealParameter(new Double[]{0.1543038});
        NodeMath nodeMath = new NodeMath();
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "shrinkage", true, "oneRateOnly", true);

        // shrinkage method
        nodeMath.estimateCorrelations(traitValuesRM);
        nodeMath.populateShrinkageEstimation(0.25925925925926);
        nodeMath.populateTraitRateMatrix();
        nodeMath.populateInverseTraitRateMatrix();
        nodeMath.populateTransformedTraitValues(traitValuesRM);
        double[] transformedTraitValues = nodeMath.getTransformedTraitValues();

        // initialize the mVector for parent of A and B
        double[] mVecD = new double [2];

        // the first node, i.e. tip A
        int tipIdx = 0;
        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(tipIdx).getLength(), nTraits, tipIdx);
        Assert.assertEquals(-0.040186380252226275, nodeMath.getAForNode(tipIdx), 0.0);
        Assert.assertEquals(-0.040186380252226275, nodeMath.getCForNode(tipIdx), 0.0);
        Assert.assertEquals(0.08037276050445255, nodeMath.getEForNode(tipIdx), 0.0);
        Assert.assertEquals(-2.11356981107731, nodeMath.getfForNode(tipIdx), EPSILON);

        PruneLikelihoodUtils.populateLmrForTipWithShrinkage(nodeMath, transformedTraitValues, nTraits, tipIdx);
        Assert.assertEquals(-0.0401863802522263,nodeMath.getLForNode(tipIdx), EPSILON);
        Assert.assertArrayEquals(new Double[] {0.20460704078553443, 0.37944671022339144}, ArrayUtils.toObject(nodeMath.getMVecForNode(tipIdx)));
        Assert.assertEquals(-3.26970682734958,nodeMath.getRForNode(tipIdx), EPSILON);

        MatrixUtilsContra.vectorAdd(mVecD, nodeMath.getMVecForNode(tipIdx), mVecD);

        // the second node, i.e. tip B
        tipIdx = 1;
        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(tipIdx).getLength(), nTraits, tipIdx);
        Assert.assertEquals(-0.0401863802522263, nodeMath.getAForNode(tipIdx), EPSILON);
        Assert.assertEquals(-0.0401863802522263, nodeMath.getCForNode(tipIdx), EPSILON);
        Assert.assertEquals(0.0803727605044526, nodeMath.getEForNode(tipIdx), EPSILON);
        Assert.assertEquals(-2.11356981107731, nodeMath.getfForNode(tipIdx), EPSILON);

        PruneLikelihoodUtils.populateLmrForTipWithShrinkage(nodeMath, transformedTraitValues, nTraits, tipIdx);
        Assert.assertEquals(-0.0401863802522263, nodeMath.getLForNode(tipIdx), EPSILON);
        Assert.assertArrayEquals(new Double[] {0.6138211223566032, 0.8401752608249581}, ArrayUtils.toObject(nodeMath.getMVecForNode(tipIdx)));
        Assert.assertEquals(-8.84887933857218, nodeMath.getRForNode(tipIdx), EPSILON);

        // mVector is added up to their parent
        MatrixUtilsContra.vectorAdd(mVecD, nodeMath.getMVecForNode(tipIdx), mVecD);

        // the third node, i.e. internal node that is the parent of A and B
        int internalIdx = 3;
        // set the L m r coming from its child nodes
        nodeMath.setLForNode(internalIdx , nodeMath.getLForNode(0) + nodeMath.getLForNode(1));
        nodeMath.setMVecForNode(internalIdx, mVecD);
        nodeMath.setRForNode(internalIdx , nodeMath.getRForNode(0) + nodeMath.getRForNode(1));
        Assert.assertEquals( -0.0803727605044526, nodeMath.getLForNode(internalIdx), EPSILON);
        Assert.assertArrayEquals(new Double[] {0.8184281631421377, 1.2196219710483496}, ArrayUtils.toObject(nodeMath.getMVecForNode(internalIdx)));
        Assert.assertEquals(-12.1185861659218, nodeMath.getRForNode(internalIdx), EPSILON);

        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(internalIdx).getLength(), nTraits, internalIdx);
        Assert.assertEquals(-0.0116480008346305, nodeMath.getAForNode(internalIdx), EPSILON);
        Assert.assertEquals(-0.0116480008346305, nodeMath.getCForNode(internalIdx), EPSILON);
        Assert.assertEquals(0.023296001669261, nodeMath.getEForNode(internalIdx), EPSILON);
        Assert.assertEquals(-3.3519633864921, nodeMath.getfForNode(internalIdx), EPSILON);

        PruneLikelihoodUtils.populateLmrForInternalNodeWithShrinkage(nodeMath, nTraits, internalIdx);
        Assert.assertEquals(-0.0101735952606144, nodeMath.getLForNode(internalIdx), EPSILON);
        Assert.assertArrayEquals(new Double[] {0.10359675130524983, 0.15437991959615796}, ArrayUtils.toObject(nodeMath.getMVecForNode(internalIdx)));
        Assert.assertEquals(-8.32455362290912, nodeMath.getRForNode(internalIdx), EPSILON);
    }

    /*
     * (4) Populate a real matrix for trait values considering the population variance
     *     Population variance is estimated by a given population, either a sample or the data itself
     *
     */
    @Test
    public void testPopulateEstimatedVarianceMatrix () {
        // trait values
        nTraits = 2;
        RealMatrix traitValuesRM = new Array2DRowRealMatrix(new double[][]{
                {-2.62762948691895, -1.56292164859448},
                {-1.50846427625826, -1.59482814741543},
                {-0.226074849617958, -2.11000367246907}
        });

        // shrinkage parameter for the population variance
        double lambda = 0.574732079259225;
        traitValuesRM = PruneLikelihoodUtils.populateTraitValueMatrixEstimatedPopulationVariance(traitValuesRM, traitValuesRM, nTraits, lambda);

        // the first species
        Assert.assertArrayEquals(new Double[] { -2.5567665115854, -2.2507927602031303 }, ArrayUtils.toObject(traitValuesRM.getRow(0)));
        // the second species
        Assert.assertArrayEquals(new Double[] { -1.4677834012215858, -2.296741907183215 }, ArrayUtils.toObject(traitValuesRM.getRow(1)));
        // the third species
        Assert.assertArrayEquals(new Double[] { -0.21997797158710666, -3.0386558368209258 }, ArrayUtils.toObject(traitValuesRM.getRow(2)));
    }

    /*
     * (5) Populate a real matrix for trait values considering the population variance
     *     Population variance is taken from input
     *
     */
    @Test
    public void testPopulateGivenPopulationVariance () {
        nTraits = 2;
        RealMatrix traitValuesRM = new Array2DRowRealMatrix(new double[][]{
                {-2.62762948691895, -1.56292164859448},
                {-1.50846427625826, -1.59482814741543},
                {-0.226074849617958, -2.11000367246907}
        });

        // the shared population variance by species
        double popVar = 0.3;
        traitValuesRM = PruneLikelihoodUtils.populateTraitValueMatrixGivenPopulationVariance(traitValuesRM, popVar, nTraits);

        // the first species
        Assert.assertArrayEquals(new Double[] { -4.7973731425041155, -2.853491475161197 }, ArrayUtils.toObject(traitValuesRM.getRow(0)));
        // the second species
        Assert.assertArrayEquals(new Double[] { -2.754066370991179, -2.911744505612018 }, ArrayUtils.toObject(traitValuesRM.getRow(1)));
        // the third species
        Assert.assertArrayEquals(new Double[] { -0.41275431606781265, -3.852322026100173 }, ArrayUtils.toObject(traitValuesRM.getRow(2)));
    }
}
