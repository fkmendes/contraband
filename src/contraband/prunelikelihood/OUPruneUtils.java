package contraband.prunelikelihood;

import beast.evolution.tree.Node;
import contraband.math.MatrixUtilsContra;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import java.util.Arrays;
import java.util.List;

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
        return eMat.preMultiply(omegaVec).mapMultiply(-1);
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
        return eMat.preMultiply(traitsValues).add(dVec);
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
        return cMat.subtract(eAPlusLInv.multiply(eMat.transpose().scalarMultiply(0.25)));
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
        return (identity.subtract(phiMat)).preMultiply(thetaVec);
    }

    public static void getOmegaVec(double[] thetaVec, double[] phiMat, double[] identity, double[] minusPhiMat, int nTraits, double[] omega) {
        MatrixUtilsContra.matrixSubtract(identity, phiMat, nTraits, minusPhiMat);
        MatrixUtilsContra.matrixPreMultiply(thetaVec, minusPhiMat, nTraits, nTraits, omega);
    }

    public static RealMatrix getPhiRM(Node aNode, RealMatrix phiMat) {
        double nodeBranchLength = aNode.getLength();
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

    // (4) matrix operations
    /*
     * This method calculates variance-covariance matrix
     * at an input internal node
     * under OU model
     * simgaRM is the evolutionary rate matrix
     * sigmaeRM is used for measurement error matrix of observed values at tips
     */
    public static RealMatrix getOUVarianceRM (Node node, RealMatrix sigmaRM, RealMatrix sigmaeRM,
                                              RealMatrix pMat, RealMatrix inverseP, EigenDecomposition decompositionH, int nTraits) {
        double nodeBranchLength = node.getLength();
        // Sigma
        RealMatrix variance = sigmaRM.multiply(sigmaRM.transpose());

        // P_1SigmaP_t = inverseP * Sigma * t(inverseP)
        variance = inverseP.multiply(variance).multiply(inverseP.transpose());

        // fLambda_ij(t) * P_1SigmaP_t
        double lambda;
        double f;
        for (int i = 0; i < nTraits; i++) {
            for (int j = 0; j < nTraits; j++) {
                if (decompositionH.getRealEigenvalue(i) == 0.0 && decompositionH.getRealEigenvalue(j) == 0.0) {
                    //fLambdaP_1SigmaP_t[i][j] = nodeHeight * P_1SigmaP_t.getEntry(i,j);
                    variance.setEntry(i, j, nodeBranchLength * variance.getEntry(i, j));
                } else {
                    lambda = decompositionH.getRealEigenvalue(i) + decompositionH.getRealEigenvalue(j);
                    f = (1 - Math.exp(-lambda * nodeBranchLength)) / lambda;
                    //fLambdaP_1SigmaP_t[i][j] = f * P_1SigmaP_t.getEntry(i,j);
                    variance.setEntry(i, j, f * variance.getEntry(i, j));
                }
            }
        }

        // variance matrix = P * (fLambda * P_1SigmaP_t) * t(P)
        variance = pMat.multiply(variance).multiply(pMat.transpose());

        if (node.isLeaf() && sigmaeRM!=null) {
            variance = variance.add(sigmaeRM.multiply(sigmaeRM.transpose()));
        }

        return variance;
    }

    public static RealMatrix getInverseVarianceRMForOU (Node node, double[] vcvMatDetArr,
                                                        RealMatrix sigmaRM, RealMatrix sigmaeRM,
                                                        RealMatrix pMat, RealMatrix inverseP, EigenDecomposition decompositionH, int nTraits) {

        // variance-covariance matrix
        RealMatrix varianceRM = getOUVarianceRM(node, sigmaRM, sigmaeRM, pMat, inverseP, decompositionH, nTraits);


        // evaluate singularity
        double[] singularValues = new SingularValueDecomposition(varianceRM).getSingularValues();
        double min = Arrays.stream(singularValues).min().getAsDouble();
        double max = Arrays.stream(singularValues).max().getAsDouble();

        EigenDecomposition decomposition = new EigenDecomposition(varianceRM);
        double[] eValues = decomposition.getRealEigenvalues();

        for (double ei : eValues) {
            if (ei < 1.0E-5) {
                SINGULAR = true;
            }
        }

        if ((min / max) < 1.0E-6) {
            SINGULAR = true;
        }


        // inverse of V and determinant of V
        // V <- V.inverse()
        try {
            LUDecomposition VMatLUD = new LUDecomposition(varianceRM);
            varianceRM = VMatLUD.getSolver().getInverse();
            vcvMatDetArr[node.getNr()] = VMatLUD.getDeterminant();
        } catch (SingularMatrixException e) {
            SINGULAR = true;
        }

        if (vcvMatDetArr[node.getNr()] == 0.0) {
            SINGULAR = true;
        }

        return varianceRM;
    }

    /*
     * This method calculates the inverse matrix of aMat + lMat, i.e. invAPlusLRM
     * and returns the product of eMat and invAplusLRM.
     * It also sets the determinant of -2 * aPlusL in log space.
     */
    public static RealMatrix getInvAPlusLRM (Node aNode, double[] negativeTwoAplusLDetArr, RealMatrix aMat, RealMatrix lMat) {
        RealMatrix AplusL = aMat.add(lMat);

        // ensure symmetry: AplusL <- 0.5 * (AplusL + AplusL.transpose)
        AplusL = (AplusL.add(AplusL.transpose())).scalarMultiply(0.5);

        try {
            // log(det(-2*AplusL))
            negativeTwoAplusLDetArr[aNode.getNr()] = Math.log(new LUDecomposition(AplusL.scalarMultiply(-2)).getDeterminant());

            // inverse of AplusL
            // AplusL <- AplusL.inverse()
            AplusL = new LUDecomposition(AplusL).getSolver().getInverse();
        }
        catch (SingularMatrixException e) {
            SINGULAR = true;
        }

        if (negativeTwoAplusLDetArr[aNode.getNr()] == 0.0) {
            SINGULAR = true;
        }

        // AplusL.inverse
        return AplusL;
    }

    // (5) PCM pruning algorithm
    public static void pruneOUPCM (Node node, int nTraits, double[] traitValuesArr,
                                   List<RealMatrix> lMatList, List<RealVector> mVecList, double[] rArr,
                                   RealMatrix sigmaRM, RealMatrix sigmaeRM,
                                   RealVector thetaVec, RealMatrix alphaMat,
                                   RealMatrix pMat, RealMatrix inverseP, EigenDecomposition decompositionH, RealMatrix identity,
                                   double[] vcvMatDetArr, double[] negativeTwoAplusLDetArr) {

        int thisNodeIdx = node.getNr();
        RealMatrix thisNodeLMat = new Array2DRowRealMatrix(new double [nTraits][nTraits]);
        RealVector thisNodeMVec = new ArrayRealVector(new double [nTraits]);
        double thisNodeR = 0;

        List<Node> children = node.getChildren();

        for (Node child : children) {
            int childIdx = child.getNr();

            // For OU, variance matrix, Phi and Omega need to be calculated for this node.
            RealMatrix phiRM = getPhiRM(child, alphaMat);

            RealVector omegaVec = getOmegaVec(thetaVec, phiRM, identity);

            // inverse of variance-covariance of this node
            RealMatrix invVCVMat = getInverseVarianceRMForOU(child, vcvMatDetArr,
                    sigmaRM, sigmaeRM,
                    pMat, inverseP, decompositionH, nTraits);

            double varianceRMDet = vcvMatDetArr[childIdx];

            // For OU
            // matrix A, C, E
            // vector b, d
            // and value f
            // need to be calculated
            // for this node
            RealMatrix aMat = getAMatForOU(invVCVMat);
            RealMatrix eMat = getEMatForOU(phiRM, invVCVMat);
            RealMatrix cMat = getCMatForOU(phiRM, eMat);
            double f = getFforOU(omegaVec, invVCVMat, varianceRMDet, nTraits);
            RealVector bVec = getBVecForOU(invVCVMat, omegaVec);
            RealVector dVec = getDVecForOU(eMat, omegaVec);

            if (child.isLeaf()) {
                // vector of trait values at this tip
                //RealVector traitsVec = traitsValuesList.get(childIdx);
                double[] traitsArr = new double [nTraits];
                MatrixUtilsContra.getMatrixRow(traitValuesArr, childIdx, nTraits, traitsArr);
                RealVector traitsVec = new ArrayRealVector(traitsArr);

                // set the L matrix
                thisNodeLMat = thisNodeLMat.add(getLMatForOULeaf(cMat));

                // set r value
                thisNodeR += getRForOULeaf(aMat, traitsVec, bVec, f);

                // set m vector
                thisNodeMVec = thisNodeMVec.add(getMVecForOULeaf(eMat, traitsVec, dVec));
            } else {

                pruneOUPCM(child, nTraits, traitValuesArr,
                        lMatList, mVecList, rArr,
                        sigmaRM, sigmaeRM,
                        thetaVec, alphaMat, pMat,
                        inverseP, decompositionH,
                        identity, vcvMatDetArr, negativeTwoAplusLDetArr);

                // (aMat + lMat).inverse
                RealMatrix aPlusLInv = getInvAPlusLRM(child, negativeTwoAplusLDetArr, aMat, lMatList.get(childIdx));

                // determinant of -2 * (aMat + lMat)
                double logDetVNode = negativeTwoAplusLDetArr[childIdx];

                // set r value
                thisNodeR += getRForOUIntNode(bVec, mVecList.get(childIdx), aPlusLInv, f, rArr[childIdx], nTraits, logDetVNode);

                RealMatrix eAPlusLInv = eMat.multiply(aPlusLInv);
                // set m vector
                thisNodeMVec = thisNodeMVec.add(getMVecForOUIntNode(eAPlusLInv, bVec, mVecList.get(childIdx), dVec));

                // set L matrix
                thisNodeLMat = thisNodeLMat.add(getLMatForOUIntNode(cMat, eMat, eAPlusLInv));
            }
        }
        lMatList.set(thisNodeIdx, thisNodeLMat);
        mVecList.set(thisNodeIdx, thisNodeMVec);
        rArr[thisNodeIdx] = thisNodeR;
    }

    // populate matrix
    public static void populateSigmaMatrix(RealMatrix rm, double[] values) {
        int k = 0;
        for (int i = 0; i < rm.getColumnDimension(); i++) {
            for (int j = i; j < rm.getColumnDimension(); j++) {
                rm.setEntry(i, j, values[k]);
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

    public static void populateRealVector(RealVector vec, int j, double[] values) {
        int nTraits = vec.getDimension();
        for (int i = 0; i < nTraits; i ++) {
            vec.setEntry(i, values[j*nTraits + i]);
        }
    }
}
