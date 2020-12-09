package contraband.prunelikelihood;

import beast.evolution.tree.Node;
import contraband.math.MatrixUtilsContra;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;


public class OUPruneUtils {
    private static double LOGTWOPI = FastMath.log(2 * Math.PI);
    static boolean SINGULAR;
    // (1) AbCdEf
    /*
     * For both BM and OU
     * aMat = -0.5 * V.inverse
     */
    public static RealMatrix getAMatForOU (RealMatrix invVarianceRM) {
        return invVarianceRM.scalarMultiply(-0.5);
    }

    public static void getAMatForOU (double[] invVarianceRM, double[] aMat, int nTraits) {
        MatrixUtilsContra.matrixScalarMultiply(invVarianceRM, -0.5, nTraits, aMat);
    }

    /*
     * For OU
     * bVec = V.inverse * omegaVec
     */
    public static RealVector getBVecForOU (RealMatrix invVarainceRM, RealVector omegaVec) {
        return invVarainceRM.preMultiply(omegaVec);
    }

    public static void getBVecForOU (double[] invVarainceRM, double[] omegaVec, int nCol, int nRow,double[] bVec) {
        MatrixUtilsContra.matrixPreMultiply(omegaVec, invVarainceRM, nCol, nRow, bVec);
    }

    /*
     * For OU
     * cMat = -0.5 * Phi.transpose * V_ * Phi
     */
    public static RealMatrix getCMatForOU (RealMatrix phiMat, RealMatrix eMat) {
        return eMat.multiply(phiMat).scalarMultiply(-0.5);
    }
    public static void getCMatForOU (double[] phiMat, double[] eMat, int nCol, int nRow, double[] cMat) {
       MatrixUtilsContra.matrixMultiplyScalar(eMat, phiMat, -0.5, nCol, nRow, cMat);
    }

    /*
     * For OU
     * dVec = - eMat * omegaVec
     */
    public static RealVector getDVecForOU (RealMatrix eMat, RealVector omegaVec) {
        RealVector a = eMat.preMultiply(omegaVec);
        RealVector b = a.mapMultiply(-1);
        RealMatrix c = eMat.scalarMultiply(-1);
        RealVector d = c.preMultiply(omegaVec);
        RealVector e = eMat.transpose().preMultiply(omegaVec).mapMultiply(-1);
        return eMat.transpose().preMultiply(omegaVec).mapMultiply(-1);
    }

    public static void getDVecForOU (double[] eMat, double[] omegaVec, int nCol, int nRow, double[] dVec) {
        MatrixUtilsContra.matrixPreMapMultiply(omegaVec, eMat, -1, nCol, nRow, dVec);
    }

    /*
     * For OU
     * eMat = Phi.transpose * V_
     */
    public static RealMatrix getEMatForOU (RealMatrix phiMat, RealMatrix invVarianceRM) {
        return phiMat.transpose().multiply(invVarianceRM);
    }
    public static void getEMatForOU (double[] phiMat, double[] invVarianceRM, int nRow, int nCol, double[] eMat) {
        MatrixUtilsContra.matrixTransMultiply(phiMat, invVarianceRM, nRow, nCol, eMat);
    }

    /*
     * For OU
     * f = -0.5 * (t(omega[ki,i]) %*% V_1[ki,ki,i] %*% omega[ki,i] + sum(ki)*log(2*pi) + log(det(as.matrix(V[ki,ki,i]))))
     * TO DO : sum(Ki), temporarily = nTraits
     */
    public static double getFforOU (RealVector omegaVec, RealMatrix varianceRM, double varianceRMDet, int nTraits) {
        return -0.5 * (varianceRM.preMultiply(omegaVec).dotProduct(omegaVec) + nTraits * LOGTWOPI +  Math.log(varianceRMDet));
    }

    public static double getFforOU (double[] omegaVec, double[] varianceRM, double varianceRMDet, int nTraits) {
        return -0.5 * (MatrixUtilsContra.tVecDotMatrixDotVec(omegaVec, varianceRM, nTraits) + nTraits * LOGTWOPI +  Math.log(varianceRMDet));
    }

    // (2) Lmr for leaf node
    /*
     * For both BM and OU
     * lMat of a lead node
     * L = 0.5 * (C + C.transpose)
     */
    public static RealMatrix getLMatForOULeaf (RealMatrix cMat){
        return  cMat.add(cMat.transpose()).scalarMultiply(0.5);
    }
    public static void getLMatForOULeaf (double[] cMat, int nTraits, double[] lMat) {
        MatrixUtilsContra.matrixTransAddScalar(cMat, 0.5, nTraits, lMat);
    }

    /*
     * For OU
     * m vector of a leaf node
     * m = d + (E * X)
     *
     */
    public static RealVector getMVecForOULeaf (RealMatrix eMat, RealVector traitsValues, RealVector dVec) {
        return eMat.transpose().preMultiply(traitsValues).add(dVec);
    }
    public static void getMVecForOULeaf (double[] eMat, double[] traitsValues, double[] dVec, int nCol, int nRow, double[] lMat) {
        MatrixUtilsContra.matrixPreMultiplyAddVector(traitsValues, eMat, dVec, nCol, nRow, lMat);
    }

    /*
     * For OU
     * r value of a leaf node
     * r = (X.transpose * A * X) + (X.transpose * b) + f
     *
     */
    public static double getRForOULeaf (RealMatrix aMat, RealVector traitsValues, RealVector bVec, double f) {
        return aMat.preMultiply(traitsValues).dotProduct(traitsValues) + traitsValues.dotProduct(bVec) + f;
    }
    public static double getRForOULeaf (double[] aMat, double[] traitsValues, double[] bVec, double f, int nTraits) {
        return MatrixUtilsContra.tVecDotMatrixDotVec(traitsValues, aMat, nTraits) + MatrixUtilsContra.vectorDotMultiply(traitsValues, bVec) + f;
    }

    // (3) Lmr for internal node
    /*
     * For both BM and OU
     * L =
     * C - 0.25 * E * (A + L).inverse * E.transpose
     */
    public static RealMatrix getLMatForOUIntNode (RealMatrix cMat, RealMatrix eMat, RealMatrix eAPlusLInv) {
        return cMat.subtract(eAPlusLInv.multiply(eMat.transpose()).scalarMultiply(0.25));
    }
    public static void getLMatForOUIntNode (double[] cMat, double[] eMat, double[] eMatTransScalar, double[] cToSubtract,double[] eAPlusLInv, int nTraits, double[] lMat) {
        MatrixUtilsContra.matrixTransScalar(eMat, 0.25, nTraits, eMatTransScalar);
        MatrixUtilsContra.matrixMultiply(eAPlusLInv, eMatTransScalar, nTraits, nTraits, cToSubtract);
        MatrixUtilsContra.matrixSubtract(cMat, cToSubtract, nTraits, lMat);
    }

    /*
     * For OU
     * m vector of an internal node
     * m = d - 0.5 * E * (A + L).inverse * （bVec + mVec)
     *
     */
    public static RealVector getMVecForOUIntNode (RealMatrix eAPlusLInv, RealVector bVec, RealVector mVec, RealVector dVec) {
        // m <- b + m
        mVec = bVec.add(mVec);

        return dVec.subtract(eAPlusLInv.transpose().preMultiply(mVec).mapMultiply(0.5));
    }

    public static void getMVecForOUIntNode (double[] eAPlusLInv, double[] bVec, double[] mCVec, double[] dToSubtract, int nTraits, double[] dVec, double[] mVec) {
        // m <- b + m
        MatrixUtilsContra.vectorAdd(bVec, mCVec, mCVec);
        MatrixUtilsContra.matrixTransPreMapMultiply(mCVec, eAPlusLInv,0.5, nTraits, nTraits, dToSubtract);
        MatrixUtilsContra.vectorSubtract(dVec, dToSubtract, mVec);
    }

    /*
     * For OU
     * r value of an internal node
     * r =
     * f + r + 0.5 * nTraits * log(2*PI)
     * - 0.5 * log(det(-2 * (A + L)))
     * - 0.25 * (b + m).transpose * (A + L).inverse * (b + m)
     *
     */
    public static double getRForOUIntNode (RealVector bVec, RealVector mVec, RealMatrix invAPlusLRM, double f, double r, int nTraits, double logDetVNode) {
        // m <- b + m
        mVec = bVec.add(mVec);

        return r + f + 0.5 * nTraits* LOGTWOPI
                - 0.5 * logDetVNode  - 0.25 * mVec.dotProduct(invAPlusLRM.preMultiply(mVec));
    }

    public static double getRForOUIntNode (double[] bVec, double[] mVec, double[] invAPlusLRM, double f, double r, int nTraits, double logDetVNode) {
        // m <- b + m
        MatrixUtilsContra.vectorAdd(bVec, mVec, mVec);

        return r + f + 0.5 * nTraits* LOGTWOPI
                - 0.5 * logDetVNode  - 0.25 * MatrixUtilsContra.tVecDotMatrixDotVec(mVec, invAPlusLRM, nTraits);
    }

    // (3) populate omega and phi
    /*
     * For OU, this method calculates omegaVec at a node.
     * I is an identity matrix.
     *
     * omega <- (I - e_Ht) %*% Theta
     */
    public static RealVector getOmegaVec(RealVector thetaVec, RealMatrix phiMat, RealMatrix identity) {
        return identity.subtract(phiMat).transpose().preMultiply(thetaVec);
    }

    public static void getOmegaVec(double[] thetaVec, double[] phiMat, double[] identity, double[] minusPhiMat, int nTraits, double[] omega) {
        MatrixUtilsContra.matrixSubtract(identity, phiMat, nTraits, minusPhiMat);
        MatrixUtilsContra.matrixPreMultiply(thetaVec, minusPhiMat, nTraits, nTraits, omega);
    }

    public static RealMatrix getPhiRM(Node aNode, RealMatrix phiMat) {
        double nodeBranchLength = aNode.getLength();
        if(nodeBranchLength == 0.0) {
            nodeBranchLength = aNode.getParent().getLength();
        }
        return exp(phiMat.scalarMultiply(-nodeBranchLength));
    }

    public static void getPhiMatArr(Node aNode, RealMatrix phiMat, int nTraits, double[] phi) {
        double nodeBranchLength = aNode.getLength();
        RealMatrix phiRM = exp(phiMat.scalarMultiply(-nodeBranchLength));
        MatrixUtilsContra.populateMatrixArray(phiRM, nTraits, nTraits, phi);
    }

    /*
     * The following two methods calculate matrix exponential of a RealMatrix,
     * by using a method in jblas that applies to DoubleMatrix.
     */
    public static RealMatrix exp(RealMatrix matrix) {
        return new Array2DRowRealMatrix(exp(matrix.getData()));
    }

    public static double[][] exp(double[][] matrix) {
        return MatrixFunctions.expm(new DoubleMatrix(matrix)).toArray2();
    }

    // populate Sigma matrix
    // this method populates an upper diagonal matrix with diagonal values and off-diagonal values
    // i.e. variance and covariance
    public static void populateUpperSigmaMatrix(RealMatrix rm, double[] diagValues, double[] offDiagValues, int nTraits) {
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            rm.setEntry(i, i, diagValues[i]);
            for (int j = i+1; j < nTraits; j++) {
                rm.setEntry(i, j, offDiagValues[k]);
                rm.setEntry(j, i, 0);
                k++;
            }
        }
    }

    // this method populates a complete matrix with evolutionary rates and correlations
    public static void populateSigmaMatrix(RealMatrix rm,  double[] diagValues, double[] offDiagValues, int nTraits) {
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            rm.setEntry(i, i, diagValues[i]);
            for (int j = i + 1; j < nTraits; j++) {
                double value = FastMath.sqrt(diagValues[i]) * FastMath.sqrt(diagValues[j]) * offDiagValues[k];
                rm.setEntry(i, j, value);
                rm.setEntry(j, i, value);
                k++;
            }
        }
    }

    // this method populates a upper diagonal matrix with evolutionary rates and correlations
    public static void populateUpperRhoSigmaMatrix(RealMatrix rm,  double[] diagValues, double[] offDiagValues, int nTraits) {
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            rm.setEntry(i, i, diagValues[i]);
            for (int j = i + 1; j < nTraits; j++) {
                double value = FastMath.sqrt(diagValues[i]) * FastMath.sqrt(diagValues[j]) * offDiagValues[k];
                rm.setEntry(i, j, value);
                rm.setEntry(j, i, 0);
                k++;
            }
        }
    }

    // this method populates a complete diagonal matrix with diagonal values and off-diagonal values
    // i.e. variance and covariance
    public static void populateDirectSigmaMatrix(RealMatrix rm,  double[] diagValues, double[] offDiagValues, int nTraits) {
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            rm.setEntry(i, i, diagValues[i]);
            for (int j = i + 1; j < nTraits; j++) {
                rm.setEntry(i, j, offDiagValues[k]);
                rm.setEntry(j, i, offDiagValues[k]);
                k++;
            }
        }
    }


    public static void populateAlphaMatrix(RealMatrix rm, double[] values) {
        int k = 0;
        for (int i = 0; i < rm.getColumnDimension(); i++) {
            for (int j = 0; j < rm.getColumnDimension(); j++) {
                rm.setEntry(i, j, values[k]);
                k++;
            }
        }
    }

    public static void populateRealVector(RealVector vec, int nTraits, double[] values) {
        for (int i = 0; i < nTraits; i ++) {
            vec.setEntry(i, values[i]);
        }
    }
}
