package contraband.utils;

import contraband.math.GeneralNodeMath;
import contraband.math.LUDecompositionForArray;
import contraband.math.MatrixUtilsContra;
import net.jsign.bouncycastle.pqc.math.linearalgebra.Matrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;

public class MorphologyLikelihoodUtils {

    private static double LOGTWOPI = FastMath.log(2 * Math.PI);


    public static void populateACEf(GeneralNodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {

        double variance = 1.0 / branchLength;

        double a = -0.5 * variance;

        // set the values in nodeMath
        nodeMath.setAForNode(nodeIdx, a);
        nodeMath.setEForNode(nodeIdx, variance);
        nodeMath.setCForNode(nodeIdx, a);
        nodeMath.setfForNode(nodeIdx, -0.5 * (nTraits * LOGTWOPI + nTraits * Math.log(branchLength) + nodeMath.getTraitRateMatrixDeterminant()));
    }

    public static void populateACEfMatrix(GeneralNodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {

        // variance = r * t * sigma + sigmaE
        nodeMath.populateVarianceMatrix(branchLength);

        // check nearly singular or not
        nodeMath.checkNearlySingularMatrix();

        // LUDecomposition -> inverse matrix and determinant
        nodeMath.operateOnVarianceMatrix();


        // calculate A, C, E, f
        //RealMatrix aMat = PCMUtils.getAMatForBM(invVCVMat); invVarianceRM.scalarMultiply(-0.5);
        //RealMatrix eMat = PCMUtils.getEMatForBM(invVCVMat); invVarianceRM
        //RealMatrix cMat = PCMUtils.getCMatForBM(eMat); eMat.scalarMultiply(-0.5)
        //double f = PCMUtils.getFforBM(varianceRMDet, nTraits); -0.5 * (nTraits * LOGTWOPI + Math.log(varianceRMDet))
        double[] res = new double[nTraits * nTraits];
        MatrixUtilsContra.vectorMapMultiply(nodeMath.getInvVarianceMatrix(), -0.5, res);
        nodeMath.setfForNode(nodeIdx, -0.5 * (nTraits * LOGTWOPI + nodeMath.getVarianceMatrixDet()));

        // set the values in nodeMath
        //nodeMath.setAForNode(nodeIdx, a);
        //nodeMath.setEForNode(nodeIdx, variance);
        //nodeMath.setCForNode(nodeIdx, a);
        //nodeMath.setfForNode(nodeIdx, -0.5 * (nTraits * LOGTWOPI + nTraits * Math.log(branchLength) + nodeMath.getTraitRateMatrixDeterminant()));
    }

    /*
     * This method calculates L m r under BM model.
     * L and r are single double values, m is a vector.
     *
     * Has side effect (sets state of nodeMath)
     */
    public static void populateLmrForTip(GeneralNodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        // vector of trait values at this tip
        // m_i <- nodeMath.traitsVec
       nodeMath.setTraitsVecForTip(traitValuesArr, nodeIdx);

        // node expectation at the tips is equal to its trait values
        nodeMath.setExpectationForTip(nodeIdx);

         // set the L matrix
        // L_tips = C_i)
        nodeMath.setLForNode(nodeIdx, new double[] {nodeMath.getCForNode(nodeIdx)});

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
    public static void populateLmrForTipTransform(GeneralNodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        // vector of trait values at this tip
        // m_i <- nodeMath.traitsVec
        nodeMath.setTraitsVecForTip(traitValuesArr, nodeIdx);

        // node expectation at the tips is equal to its trait values
        nodeMath.setExpectationForTip(nodeIdx);

        // set the L matrix
        // L_tips = C_i
        nodeMath.setLForNode(nodeIdx, new double[] {nodeMath.getCForNode(nodeIdx)});

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
    public static void populateLmrForInternalNode(GeneralNodeMath nodeMath, int nTraits, int nodeIdx) {
        // (A + L)
        double aPlusL = nodeMath.getAForNode(nodeIdx)  + nodeMath.getLForNode(nodeIdx)[0];

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
        nodeMath.setLForNode(nodeIdx,new double[] {nodeMath.getCForNode(nodeIdx)  - 0.25 * eAPlusLInv * nodeMath.getEForNode(nodeIdx)});
    }

    public static void populateLmrForInternalNodeTransform(GeneralNodeMath nodeMath, int nTraits, int nodeIdx) {
        // (A + L)
        double aPlusL = nodeMath.getAForNode(nodeIdx)  + nodeMath.getLForNode(nodeIdx)[0];

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
        nodeMath.setLForNode(nodeIdx,new double[]{nodeMath.getCForNode(nodeIdx)  - 0.25 * eAPlusLInv * nodeMath.getEForNode(nodeIdx)});
    }
}
