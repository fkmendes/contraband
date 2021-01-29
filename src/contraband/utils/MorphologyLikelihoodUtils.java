package contraband.utils;

import contraband.math.GeneralNodeMath;
import contraband.math.MatrixUtilsContra;
import contraband.prunelikelihood.MorphologicalData;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.util.FastMath;
import outercore.parameter.KeyRealParameter;

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

    public static void populateACEfMatrixForTips(GeneralNodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {
        // variance = r * t * sigma + sigmaE
        nodeMath.populateVarianceMatrixForTips(branchLength);

        populateACEFMatrix(nodeMath, nTraits, nodeIdx);
    }

    public static void populateACEfMatrixForIntNodes(GeneralNodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {
        // variance = r * t * sigma + sigmaE
        nodeMath.populateVarianceMatrixForIntNodes(branchLength);

        populateACEFMatrix(nodeMath, nTraits, nodeIdx);
    }

    public static void populateACEFMatrix(GeneralNodeMath nodeMath, int nTraits, int nodeIdx) {
        // check nearly singular or not
        nodeMath.checkNearlySingularMatrix();

        // LUDecomposition -> inverse matrix and determinant
        nodeMath.operateOnVarianceMatrix();


        // calculate A, C, E, f
        // A = -0.5 * V.inverse
        // C = -0.5 * V.inverse
        // E = -V.inverse
        // f = -0.5 * (nTraits * LOGTWOPI + log(V.determinant))
        double[] res = new double[nTraits * nTraits];
        MatrixUtilsContra.vectorMapMultiply(nodeMath.getInvVarianceMatrix(), -0.5, res);

        nodeMath.setAMatForNode(nodeIdx, res);

        nodeMath.setCMatForNode(nodeIdx, res);

        nodeMath.setEMatForNode(nodeIdx, nodeMath.getInvVarianceMatrix());

        double f = -0.5 * (nTraits * LOGTWOPI + nodeMath.getVarianceMatrixDet());
        nodeMath.setfForNode(nodeIdx, f);
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
        nodeMath.setLForNode(nodeIdx, nodeMath.getCMatForNode(nodeIdx));

        // set r value
        // r_tip = A_i * (m_i.transpose * invTraitRateMatrix * m_i) + f_i
        nodeMath.setRForNode(nodeIdx, nodeMath.getAMatForNode(nodeIdx)[0] *
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
                    nodeMath.getEMatForNode(nodeIdx)[0],
                    nodeMath.getTempVec());

        nodeMath.setMVecForNode(nodeIdx, nodeMath.getTempVec());
    }

    public static void populateLmrMatrixForTip(GeneralNodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        // vector of trait values at this tip
        // m_i <- nodeMath.traitsVec
        nodeMath.setTraitsVecForTip(traitValuesArr, nodeIdx);

        // set the L matrix
        // L_tips = C_i)
        nodeMath.setLForNode(nodeIdx, nodeMath.getCMatForNode(nodeIdx));

        // set r value
        // r_tip = (m_i.transpose * A * m_i) + f_i
        double r = MatrixUtilsContra.tVecDotMatrixDotVec(nodeMath.getTraitsVec(), nodeMath.getAMatForNode(nodeIdx), nTraits) + nodeMath.getfForNode(nodeIdx);
        nodeMath.setRForNode(nodeIdx, r);

        // calculate m vector
        // m_tip = E * m_i
        // TO DO: check  E %*% m_i and m_i * E
        MatrixUtilsContra.matrixPreMultiply(nodeMath.getTraitsVec(),
                nodeMath.getEMatForNode(nodeIdx),
                nTraits,
                nTraits,
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
        nodeMath.setLForNode(nodeIdx, nodeMath.getCMatForNode(nodeIdx));

        // set r value
        // r_tip = A_i * (m_i.transpose * invTraitRateMatrix * m_i) + f_i
        nodeMath.setRForNode(nodeIdx, MatrixUtilsContra.vecTransScalarMultiply(
                    nodeMath.getTraitsVec(),
                    nodeMath.getAMatForNode(nodeIdx)[0],
                    nTraits) +
                    nodeMath.getfForNode(nodeIdx));

        // calculate m vector
        // m_tip = E * invTraitRate * m_i
            MatrixUtilsContra.vectorMapMultiply(nodeMath.getTraitsVec(),
                    nodeMath.getEMatForNode(nodeIdx)[0],
                    nodeMath.getTempVec());

        nodeMath.setMVecForNode(nodeIdx, nodeMath.getTempVec());
    }

    /*
     * This method calculates L m r under BM model.
     * L and r are single double values, m is a vector.
     *
     * Has side effect (sets state of nodeMath)
     */
    public static void populateLmrForIntNode(GeneralNodeMath nodeMath, int nTraits, int nodeIdx) {
        // (A + L)
        double aPlusL = nodeMath.getAMatForNode(nodeIdx)[0] + nodeMath.getLMatForNode(nodeIdx)[0];

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
        double eAPlusLInv = nodeMath.getEMatForNode(nodeIdx)[0] * aPlusLInv;
        MatrixUtilsContra.vectorMapMultiply(mVecForChild, -0.5 * eAPlusLInv, nodeMath.getTempVec());
        nodeMath.setMVecForNode(nodeIdx, nodeMath.getTempVec());

        // set L value
        // L = C - 0.25 * E * (A + L).inverse * E.transpose
        nodeMath.setLForNode(nodeIdx,new double[] {nodeMath.getCMatForNode(nodeIdx)[0]  - 0.25 * eAPlusLInv * nodeMath.getEMatForNode(nodeIdx)[0]});
    }

    public static void populateLmrMatrixForIntNode(GeneralNodeMath nodeMath, int nTraits, int nodeIdx) {
        // get (A + L).inverse and |(-2) * (A + L)|
        nodeMath.operateOnAPlusLMatrixForNode(nodeIdx);

        // set r value
        // r_non-tip = f_i + r_i + (nTraits / 2) * log(2 * pi) - (1/2) * log|-2 * (A_i + L_i)| - (1/4) * m_i.transpose * (A + L).inverse * m_i
        double [] mVecForChild = nodeMath.getMVecForNode(nodeIdx);
        nodeMath.setRForNode(nodeIdx, nodeMath.getRForNode(nodeIdx) +
                nodeMath.getfForNode(nodeIdx) +
                0.5 * nTraits* LOGTWOPI - 0.5 * nodeMath.getLogDetNegativeTwoAplusL() -
                0.25 * MatrixUtilsContra.tVecDotMatrixDotVec(
                        mVecForChild,
                        nodeMath.getAPlusLInverse(),
                        nTraits));

        // set m vector
        // m_non-tip = - 0.5 *  E * (A + L).inverse * m_i
        double[] eAPlusLInv = new double[nTraits * nTraits];
        double[] eAplusLInvTrans = new double[nTraits * nTraits];
        MatrixUtilsContra.matrixMultiply(nodeMath.getEMatForNode(nodeIdx), nodeMath.getAPlusLInverse(), nTraits, nTraits, eAplusLInvTrans);
        MatrixUtilsContra.matrixTranspose(eAplusLInvTrans, nTraits, eAPlusLInv);

        double[] res = new double[nTraits];
        MatrixUtilsContra.matrixPreMultiply(mVecForChild, eAPlusLInv, nTraits, nTraits, res);
        MatrixUtilsContra.vectorMapMultiply(res, -0.5, nodeMath.getTempVec());
        nodeMath.setMVecForNode(nodeIdx, nodeMath.getTempVec());

        // set L value
        // L = C - 0.25 * E * (A + L).inverse * E.transpose
        double[] eTrans = new double [nTraits * nTraits];
        MatrixUtilsContra.matrixTranspose(nodeMath.getEMatForNode(nodeIdx), nTraits, eTrans);
        double[] eTransScaler =  new double [nTraits * nTraits];
        MatrixUtilsContra.vectorMapMultiply(eTrans, -0.25, eTransScaler);
        double[] mat = new double [nTraits * nTraits];
        MatrixUtilsContra.matrixMultiply(eAplusLInvTrans, eTransScaler, nTraits, nTraits, mat);
        double[] res1 = new double [nTraits * nTraits];
        MatrixUtilsContra.vectorAdd(nodeMath.getCMatForNode(nodeIdx), mat, res1);
        nodeMath.setLForNode(nodeIdx, res1);
    }

    public static void populateLmrForInternalNodeTransform(GeneralNodeMath nodeMath, int nTraits, int nodeIdx) {
        // (A + L)
        double aPlusL = nodeMath.getAMatForNode(nodeIdx)[0]  + nodeMath.getLMatForNode(nodeIdx)[0];

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
        double eAPlusLInv = nodeMath.getEMatForNode(nodeIdx)[0] * aPlusLInv;
        MatrixUtilsContra.vectorMapMultiply(mVecForChild, -0.5 * eAPlusLInv, nodeMath.getTempVec());
        nodeMath.setMVecForNode(nodeIdx, nodeMath.getTempVec());

        // set L value
        // L = C - 0.25 * E * (A + L).inverse * E.transpose
        nodeMath.setLForNode(nodeIdx,new double[]{nodeMath.getCMatForNode(nodeIdx)[0]  - 0.25 * eAPlusLInv * nodeMath.getEMatForNode(nodeIdx)[0]});
    }


    public static RealMatrix populateTraitMatrixForPopulationSample(KeyRealParameter traitsValues){
        int nTraits = traitsValues.getMinorDimension1();
        int nSpecies = traitsValues.getMinorDimension2();
        RealMatrix traitRM = new Array2DRowRealMatrix(new double[nSpecies][nTraits]);
        String[] keys = traitsValues.getKeys();
        for (int i = 0; i < nSpecies; i ++) {
            // get all traits values for this species
            Double[] traitForSpecies = traitsValues.getRowValues(keys[i]);
            for (int j= 0; j < nTraits; j ++) {
                traitRM.setEntry(i, j, traitForSpecies[j]);
            }
        }
        return traitRM;
    }

    public static RealMatrix standardiseContinuousTraits(RealMatrix popRM, RealMatrix contRM, int nTraits, double lambda) {
        Variance var = new Variance();
        Median median= new Median();
        double[] empVariance = new double[nTraits];
        double [] empStandardDeviation = new double[nTraits];

        for(int i = 0; i < nTraits; i++) {
            empVariance[i] = var.evaluate(popRM.getColumn(i));
        }

        double target = median.evaluate(empVariance);
        for (int k = 0; k < nTraits; k ++) {
            empStandardDeviation[k] = 1.0 / Math.sqrt(lambda * target + (1 - lambda) * empVariance[k]);
        }

        return contRM.multiply(MatrixUtils.createRealDiagonalMatrix(empStandardDeviation));
    }

    public static RealMatrix getContTraitRealMatrix (double[] traitValuesArr, int nTraits, int nSpecies) {
        RealMatrix traitRM = new Array2DRowRealMatrix(new double[nSpecies][nTraits]);
        for (int i = 0; i < nSpecies; i ++ ) {
            for (int j = 0; j < nTraits; j++) {
                traitRM.setEntry(i, j, traitValuesArr[i * nTraits + j]);
            }
        }
        return traitRM;
    }

    // This method calculates unbiased estimation of correlation matrix
    public static void populateUnbiasedRho(RealMatrix traitMat, double[] unbiasedRho){
        PearsonsCorrelation correlation = new PearsonsCorrelation(traitMat);
        RealMatrix unbiasedRhoRM = correlation.getCorrelationMatrix();

        int idx = 0;
        int dim = unbiasedRhoRM.getColumnDimension();
        for (int i = 0; i < dim; i++){
            for (int j = i + 1; j < dim; j ++) {
                unbiasedRho[idx] =  unbiasedRhoRM.getEntry(i, j);
                idx ++;
            }
        }
    }
}