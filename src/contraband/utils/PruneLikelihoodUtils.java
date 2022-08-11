package contraband.utils;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Tree;
import contraband.math.NodeMath;
import contraband.math.MatrixUtilsContra;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.util.FastMath;


public class PruneLikelihoodUtils {

    private static double LOGTWOPI = FastMath.log(2 * Math.PI);

    public static void populateTraitValuesArr(RealParameter traitValues, Tree tree, int nTraits, double[] traitValuesArr) {
        // according to node number of tips
        for (int i = 0; i < tree.getLeafNodeCount(); i ++) {
            // get all traits values for this species
            Double[] traitForSpecies = traitValues.getRowValues(tree.getNode(i).getID());
            for (int j= 0; j < nTraits; j ++) {
                // populate the traits one by one in an array
                traitValuesArr[i*nTraits + j] = traitForSpecies[j];
            }
        }
    }

    public static void  populateTraitValuesArr(RealParameter traitValues, Tree tree, NodeMath nodeMath, int nTraits, double[] traitValuesArr) {
        for (int i = 0; i < tree.getLeafNodeCount(); i ++) {
            // get all traits values for this species
            Double[] traitForSpecies = traitValues.getRowValues(tree.getNode(i).getID());
            if(traitForSpecies == null) {
                nodeMath.setSpeciesToIgnore(i);
            } else {
                for (int j = 0; j < nTraits; j++) {
                    // populate the traits one by one in an array
                    traitValuesArr[i * nTraits + j] = traitForSpecies[j];
                }
            }
        }
    }

    /*
     * Fills out traitRM in place
     */
    public static void populateTraitValuesMatrix(RealParameter traitValues, Tree tree, int nTraits, RealMatrix traitRM){
        int index = 0;
        for (int i = 0; i < tree.getLeafNodeCount(); i ++) {
            // get all traits values for this species
            Double[] traitForSpecies = traitValues.getRowValues(tree.getNode(i).getID());
            if(traitForSpecies != null) {
                for (int j = 0; j < nTraits; j++) {
                    traitRM.setEntry(index, j, traitForSpecies[j]);
                }
                index ++;
            }
        }
    }

    /*
     * This method calculates A C E f under BM model.
     * A C E and f are all scalars
     *
     * Has side effect (sets state of nodeMath)
     */
    public static void populateACEf(NodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {

        double variance = 1.0 / branchLength;

        double a = -0.5 * variance;

        // set the values in nodeMath
        nodeMath.setAForNode(nodeIdx, a);
        nodeMath.setEForNode(nodeIdx, variance);
        nodeMath.setCForNode(nodeIdx, a);
        nodeMath.setfForNode(nodeIdx, -0.5 * (nTraits * LOGTWOPI + nTraits * Math.log(branchLength) + nodeMath.getTraitRateMatrixDeterminant()));
    }

    /*
     * This method calculates L m r under BM model.
     * L and r are single double values, m is a vector.
     *
     * Has side effect (sets state of nodeMath)
     */
    public static void populateLmrForTip(NodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        // vector of trait values at this tip
        // m_i <- nodeMath.traitsVec
       nodeMath.setTraitsVecForTip(traitValuesArr, nodeIdx);

        // node expectation at the tips is equal to its trait values
        nodeMath.setExpectationForTip(nodeIdx);

         // set the L matrix
        // L_tips = C_i)
        nodeMath.setLForNode(nodeIdx, nodeMath.getCForNode(nodeIdx));

        // set r value
        // r_tip = A_i * (m_i.transpose * invTraitRateMatrix * m_i) + f_i
            nodeMath.setRForNode(nodeIdx, nodeMath.getAForNode(nodeIdx) *
                    MatrixUtilsContra.tVecDotMatrixDotVec(
                            nodeMath.getTraitsVec(),
                            nodeMath.getTraitRateMatrixInverse(),
                            nTraits) +
                    nodeMath.getfForNode(nodeIdx));

        // calculate m vector
        // m_tip = E * invTraitRate * m_i
        MatrixUtilsContra.matrixPreMultiply(nodeMath.getTraitsVec(),
                    nodeMath.getTraitRateMatrixInverse(),
                    nTraits,
                    nTraits,
                    nodeMath.getTempVec());
        MatrixUtilsContra.vectorMapMultiply(nodeMath.getTempVec(),
                    nodeMath.getEForNode(nodeIdx),
                    nodeMath.getTempVec());

        nodeMath.setMVecForNode(nodeIdx, nodeMath.getTempVec());
    }

    /*
     * This method calculates L m r under BM model.
     * L and r are single double values, m is a vector.
     *
     * Has side effect (sets state of nodeMath)
     */
    public static void populateLmrForTipWithShrinkage(NodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        // vector of trait values at this tip
        // m_i <- nodeMath.traitsVec
        nodeMath.setTraitsVecForTip(traitValuesArr, nodeIdx);

        // node expectation at the tips is equal to its trait values
        nodeMath.setExpectationForTip(nodeIdx);

        // set the L matrix
        // L_tips = C_i
        nodeMath.setLForNode(nodeIdx, nodeMath.getCForNode(nodeIdx));

        // set r value
        // r_tip = A_i * (m_i.transpose * invTraitRateMatrix * m_i) + f_i
        nodeMath.setRForNode(nodeIdx, MatrixUtilsContra.vecTransScalarMultiply(
                    nodeMath.getTraitsVec(),
                    nodeMath.getAForNode(nodeIdx),
                    nTraits) +
                    nodeMath.getfForNode(nodeIdx));

        // calculate m vector
        // m_tip = E * invTraitRate * m_i
            MatrixUtilsContra.vectorMapMultiply(nodeMath.getTraitsVec(),
                    nodeMath.getEForNode(nodeIdx),
                    nodeMath.getTempVec());

        nodeMath.setMVecForNode(nodeIdx, nodeMath.getTempVec());
    }

    /*
     * This method calculates L m r under BM model.
     * L and r are single double values, m is a vector.
     *
     * Has side effect (sets state of nodeMath)
     */
    public static void populateLmrForInternalNode(NodeMath nodeMath, int nTraits, int nodeIdx) {
        // (A + L)
        double aPlusL = nodeMath.getAForNode(nodeIdx)  + nodeMath.getLForNode(nodeIdx);

        // 1 / (A + L) = (A + L).inverse
        double aPlusLInv = 1.0 / aPlusL;

        // determinant of -2 * (A + L) = ((-2 * (A + L))^K) * detInvTraitRateMat
        double logDetVNode = 0.5 * nTraits * FastMath.log(- 2 * aPlusL) + 0.5 * FastMath.log(nodeMath.getTraitRateMatrixInverseDeterminant());

        // set r value
        // r_non-tip = f_i + r_i + (nTraits / 2) * log(2 * pi) - (1/2) * log|-2 * (A_i + L_i)| - (1/4) * m_i.transpose * (A + L).inverse * m_i
        double [] mVecForChild = nodeMath.getMVecForNode(nodeIdx);
        nodeMath.setRForNode(nodeIdx, nodeMath.getRForNode(nodeIdx) +
                    nodeMath.getfForNode(nodeIdx) +
                    0.5 * nTraits* LOGTWOPI - logDetVNode -
                    0.25 * aPlusLInv * MatrixUtilsContra.tVecDotMatrixDotVec(
                            mVecForChild,
                            nodeMath.getTraitRateMatrix(),
                            nTraits));

        // set m vector
        // m_non-tip = - 0.5 *  E * (A + L).inverse * m_i
        double eAPlusLInv = nodeMath.getEForNode(nodeIdx) * aPlusLInv;
        MatrixUtilsContra.vectorMapMultiply(mVecForChild, -0.5 * eAPlusLInv, nodeMath.getTempVec());
        nodeMath.setMVecForNode(nodeIdx, nodeMath.getTempVec());

        // set L value
        // L = C - 0.25 * E * (A + L).inverse * E.transpose
        nodeMath.setLForNode(nodeIdx,nodeMath.getCForNode(nodeIdx)  - 0.25 * eAPlusLInv * nodeMath.getEForNode(nodeIdx));
    }

    public static void populateLmrForInternalNodeWithShrinkage(NodeMath nodeMath, int nTraits, int nodeIdx) {
        // (A + L)
        double aPlusL = nodeMath.getAForNode(nodeIdx)  + nodeMath.getLForNode(nodeIdx);

        // 1 / (A + L) = (A + L).inverse
        double aPlusLInv = 1.0 / aPlusL;

        // determinant of -2 * (A + L) = ((-2 * (A + L))^K) * detInvTraitRateMat
        double logDetVNode = nTraits * FastMath.log(-2 * aPlusL) + nodeMath.getTraitRateMatrixInverseDeterminant();

        // set r value
        // r_non-tip = f_i + r_i + (nTraits / 2) * log(2 * pi) - (1/2) * log|-2 * (A_i + L_i)| - (1/4) * m_i.transpose * (A + L).inverse * m_i
        double [] mVecForChild = nodeMath.getMVecForNode(nodeIdx);
        nodeMath.setRForNode(nodeIdx, nodeMath.getRForNode(nodeIdx) +
                    nodeMath.getfForNode(nodeIdx) +
                    0.5 * nTraits* LOGTWOPI - 0.5 * logDetVNode -
                    0.25 * MatrixUtilsContra.vecTransScalarMultiply(
                            mVecForChild,
                            aPlusLInv,
                            nTraits));

        // set m vector
        // m_non-tip = - 0.5 *  E * (A + L).inverse * m_i
        double eAPlusLInv = nodeMath.getEForNode(nodeIdx) * aPlusLInv;
        MatrixUtilsContra.vectorMapMultiply(mVecForChild, -0.5 * eAPlusLInv, nodeMath.getTempVec());
        nodeMath.setMVecForNode(nodeIdx, nodeMath.getTempVec());

        // set L value
        // L = C - 0.25 * E * (A + L).inverse * E.transpose
        nodeMath.setLForNode(nodeIdx,nodeMath.getCForNode(nodeIdx)  - 0.25 * eAPlusLInv * nodeMath.getEForNode(nodeIdx));
    }

    public static RealMatrix populateTraitValueMatrixEstimatedPopulationVariance(RealMatrix traitRM, RealMatrix originRM, int nTraits, double lambda) {
        Variance var = new Variance();
        Median median= new Median();
        double[] empVariance = new double[nTraits];
        double [] empStandardDeviation = new double[nTraits];

        for(int i = 0; i < nTraits; i++) {
            empVariance[i] = var.evaluate(traitRM.getColumn(i));
        }

        double target = median.evaluate(empVariance);
        for (int k = 0; k < nTraits; k ++) {
            empStandardDeviation[k] = 1.0 / Math.sqrt(lambda * target + (1 - lambda) * empVariance[k]);
        }

        return originRM.multiply(MatrixUtils.createRealDiagonalMatrix(empStandardDeviation));
    }

    public static RealMatrix populateTraitValueMatrixGivenPopulationVariance(RealMatrix originRM, double popVar, int nTraits){
        double [] empStandardDeviation = new double[nTraits];
        for (int k = 0; k < nTraits; k ++) {
            empStandardDeviation[k] = 1.0 / Math.sqrt(popVar);
        }
        return originRM.multiply(MatrixUtils.createRealDiagonalMatrix(empStandardDeviation));
    }
}
