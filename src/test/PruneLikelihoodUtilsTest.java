package test;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.math.MatrixUtilsContra;
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
 * This class contains unit tests for PruneLikelihoodUtils,
 * which populates continuous trait values and parameters for pruning likelihood calculation.
 */

public class PruneLikelihoodUtilsTest {

    private TreeParser tree;
    private String treeStr;
    private Integer nTraits;
    private String spNames;
    private List<Double> data;
    private final KeyRealParameter traitValues = new KeyRealParameter();

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
        double[] traitValuesArrayList = new double[3 * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitValues, tree, nTraits, traitValuesArrayList);

        // BM model parameters
        // When shrinkage method is not used, correlations among traits are taken from input
        RealParameter sigmasq = new RealParameter(new Double[]{0.3145740});
        RealParameter correlation = new RealParameter(new Double[]{-0.632620487603683});
        RealParameter rootValues = new RealParameter(new Double[] {-1.31465955080609, -1.79611255566212});
        NodeMath nodeMath = new NodeMath();
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "correlation", correlation, "rootValues", rootValues,"oneRateOnly", true);

        // the pruning algorithm starts from tips to the root
        // the first node, i.e. tip A
        int tipIdx = 0;
        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(tipIdx).getLength(), nTraits, tipIdx);
        Assert.assertEquals(-0.021733632865102357, nodeMath.getAForNode(tipIdx), 0.0);
        Assert.assertEquals(-0.021733632865102357, nodeMath.getCForNode(tipIdx), 0.0);
        Assert.assertEquals(0.043467265730204714, nodeMath.getEForNode(tipIdx), 0.0);
        Assert.assertEquals(-4.973624202525401, nodeMath.getfForNode(tipIdx), 0.0);

        PruneLikelihoodUtils.populateLmrForTip(nodeMath, traitValuesArrayList, nTraits, tipIdx);
        Assert.assertEquals(-0.021733632865102357,nodeMath.getLForNode(tipIdx), 0.0);
        Assert.assertArrayEquals(new Double[] {0.0, 0.0}, ArrayUtils.toObject(nodeMath.getMVecForNode(tipIdx)));
        Assert.assertEquals(-4.973624202525401,nodeMath.getRForNode(tipIdx), 0.0);

        // the second node, i.e. tip B
        tipIdx = 1;
        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(tipIdx).getLength(), nTraits, tipIdx);
        Assert.assertEquals(-0.021733632865102357, nodeMath.getAForNode(tipIdx), 0.0);
        Assert.assertEquals(-0.021733632865102357, nodeMath.getCForNode(tipIdx), 0.0);
        Assert.assertEquals(0.043467265730204714, nodeMath.getEForNode(tipIdx), 0.0);
        Assert.assertEquals(-4.973624202525401, nodeMath.getfForNode(tipIdx), 0.0);

        PruneLikelihoodUtils.populateLmrForTip(nodeMath, traitValuesArrayList, nTraits, tipIdx);
        Assert.assertEquals(-0.021733632865102357,nodeMath.getLForNode(tipIdx), 0.0);
        Assert.assertArrayEquals(new Double[] {0.0, 0.0}, ArrayUtils.toObject(nodeMath.getMVecForNode(tipIdx)));
        Assert.assertEquals(-4.973624202525401,nodeMath.getRForNode(tipIdx), 0.0);

        // the third node, i.e. internal node that is the parent of A and B
        int internalIdx = 3;
        nodeMath.setLForNode(internalIdx , nodeMath.getLForNode(0) + nodeMath.getLForNode(1));
        MatrixUtilsContra.vectorAdd(nodeMath.getMVecForNode(0), nodeMath.getMVecForNode(1), nodeMath.getTempVec());
        nodeMath.setMVecForNode(internalIdx , nodeMath.getTempVec());
        nodeMath.setRForNode(internalIdx , nodeMath.getRForNode(0) + nodeMath.getRForNode(0));

        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(internalIdx).getLength(), nTraits, internalIdx);
        Assert.assertEquals(-0.034840896606782364, nodeMath.getAForNode(internalIdx), 0.0);
        Assert.assertEquals(-0.034840896606782364, nodeMath.getCForNode(internalIdx), 0.0);
        Assert.assertEquals(0.06968179321356473, nodeMath.getEForNode(internalIdx), 0.0);
        Assert.assertEquals(-4.501693278196102, nodeMath.getfForNode(internalIdx), 0.0);

        PruneLikelihoodUtils.populateLmrForInternalNode(nodeMath, nTraits, tipIdx);
        Assert.assertEquals(-0.043467265730204714, nodeMath.getLForNode(internalIdx), 0.0);
        Assert.assertArrayEquals(new Double[] {0.0, 0.0}, ArrayUtils.toObject(nodeMath.getMVecForNode(internalIdx)));
        Assert.assertEquals(-9.947248405050802, nodeMath.getRForNode(internalIdx), 0.0);
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
        double[] traitValuesArrayList = new double[3 * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitValues, tree, nTraits, traitValuesArrayList);

        // BM model parameters
        // When shrinkage method is used, correlations among traits are estimated in NodeMath
        RealParameter sigmasq = new RealParameter(new Double[]{0.1543038});
        NodeMath nodeMath = new NodeMath();
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "shrinkage", true, "oneRateOnly", true);

        // the pruning algorithm starts from tips to the root
        // the first node, i.e. tip A
        int tipIdx = 0;
        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(tipIdx).getLength(), nTraits, tipIdx);
        Assert.assertEquals(-0.040186380252226275, nodeMath.getAForNode(tipIdx), 0.0);
        Assert.assertEquals(-0.040186380252226275, nodeMath.getCForNode(tipIdx), 0.0);
        Assert.assertEquals(0.08037276050445255, nodeMath.getEForNode(tipIdx), 0.0);
        Assert.assertEquals(-4.358957026308008, nodeMath.getfForNode(tipIdx), 0.0);

        PruneLikelihoodUtils.populateLmrForTipWithShrinkage(nodeMath, traitValuesArrayList, nTraits, tipIdx);
        Assert.assertEquals(-0.040186380252226275,nodeMath.getLForNode(tipIdx), 0.0);
        Assert.assertArrayEquals(new Double[] {0.08037276050445255, 0.1607455210089051}, ArrayUtils.toObject(nodeMath.getMVecForNode(tipIdx)));
        Assert.assertEquals(-4.55988892756914,nodeMath.getRForNode(tipIdx), 0.0);

        // the second node, i.e. tip B
        tipIdx = 1;
        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(tipIdx).getLength(), nTraits, tipIdx);
        Assert.assertEquals(-0.040186380252226275, nodeMath.getAForNode(tipIdx), 0.0);
        Assert.assertEquals(-0.040186380252226275, nodeMath.getCForNode(tipIdx), 0.0);
        Assert.assertEquals(0.08037276050445255, nodeMath.getEForNode(tipIdx), 0.0);
        Assert.assertEquals(-4.358957026308008, nodeMath.getfForNode(tipIdx), 0.0);

        PruneLikelihoodUtils.populateLmrForTipWithShrinkage(nodeMath, traitValuesArrayList, nTraits, tipIdx);
        Assert.assertEquals(-0.040186380252226275, nodeMath.getLForNode(tipIdx), 0.0);
        Assert.assertArrayEquals(new Double[] {0.24111828151335765, 0.4018638025222627}, ArrayUtils.toObject(nodeMath.getMVecForNode(tipIdx)));
        Assert.assertEquals(-5.725293954883702, nodeMath.getRForNode(tipIdx), 0.0);

        // the third node, i.e. internal node that is the parent of A and B
        int internalIdx = 3;
        nodeMath.setLForNode(internalIdx , nodeMath.getLForNode(0) + nodeMath.getLForNode(1));
        MatrixUtilsContra.vectorAdd(nodeMath.getMVecForNode(0), nodeMath.getMVecForNode(1), nodeMath.getTempVec());
        nodeMath.setMVecForNode(internalIdx , nodeMath.getTempVec());
        nodeMath.setRForNode(internalIdx , nodeMath.getRForNode(0) + nodeMath.getRForNode(0));

        PruneLikelihoodUtils.populateACEf(nodeMath, tree.getNode(internalIdx).getLength(), nTraits, internalIdx);
        Assert.assertEquals(-0.011648000834630511, nodeMath.getAForNode(internalIdx), 0.0);
        Assert.assertEquals(-0.011648000834630511, nodeMath.getCForNode(internalIdx), 0.0);
        Assert.assertEquals(0.023296001669261022, nodeMath.getEForNode(internalIdx), 0.0);
        Assert.assertEquals(-5.5973506017228045, nodeMath.getfForNode(internalIdx), 0.0);

        PruneLikelihoodUtils.populateLmrForInternalNodeWithShrinkage(nodeMath, nTraits, tipIdx);
        Assert.assertEquals(-0.08037276050445255, nodeMath.getLForNode(internalIdx), 0.0);
        Assert.assertArrayEquals(new Double[] {0.4822365630267153, 0.8037276050445255}, ArrayUtils.toObject(nodeMath.getMVecForNode(internalIdx)));
        Assert.assertEquals(-9.11977785513828, nodeMath.getRForNode(internalIdx), 0.0);
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
        double lambda = 0.4;
        traitValuesRM = PruneLikelihoodUtils.populateTraitValueMatrixEstimatedPopulationVariance(traitValuesRM, traitValuesRM, nTraits, lambda);

        // the first species
        Assert.assertArrayEquals(new Double[] { -2.42497174189288, -2.5896391144824498 }, ArrayUtils.toObject(traitValuesRM.getRow(0)));
        // the second species
        Assert.assertArrayEquals(new Double[] { -1.3921229236433845, -2.6425056912729277 }, ArrayUtils.toObject(traitValuesRM.getRow(1)));
        // the third species
        Assert.assertArrayEquals(new Double[] { -0.20863867017988777, -3.496111303366598 }, ArrayUtils.toObject(traitValuesRM.getRow(2)));
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
