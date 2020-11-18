package contraband.prunelikelihood;

import beast.evolution.tree.Node;
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

    /*
     * For OU
     * bVec = V.inverse * omegaVec
     */
    public static RealVector getBVecForOU (RealMatrix invVarainceRM, RealVector omegaVec) {
        return invVarainceRM.preMultiply(omegaVec);
    }

    /*
     * For OU
     * cMat = Phi.transpose * V_
     */
    public static RealMatrix getCMatForOU (RealMatrix phiMat, RealMatrix eMat) {
        return eMat.multiply(phiMat).scalarMultiply(-0.5);
    }

    /*
     * For OU
     * dVec = - eMat * omegaVec
     */
    public static RealVector getDVecForOU (RealMatrix eMat, RealVector omegaVec) {
        return eMat.preMultiply(omegaVec).mapMultiply(-1);
    }

    /*
     * For OU
     * eMat = Phi.transpose * V_
     */
    public static RealMatrix getEMatForOU (RealMatrix phiMat, RealMatrix invVarianceRM) {
        return phiMat.transpose().multiply(invVarianceRM);
    }

    /*
     * For OU
     * f = -0.5 * (t(omega[ki,i]) %*% V_1[ki,ki,i] %*% omega[ki,i] + sum(ki)*log(2*pi) + log(det(as.matrix(V[ki,ki,i]))))
     * TO DO : sum(Ki), temporarily = nTraits
     */
    public static double getFforOU (RealVector omegaVec, RealMatrix varianceRM, double varianceRMDet, int nTraits) {
        return -0.5 * (varianceRM.preMultiply(omegaVec).dotProduct(omegaVec) + nTraits * LOGTWOPI +  Math.log(varianceRMDet));
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

    /*
     * For OU
     * m vector of a leaf node
     * m = d + (E * X)
     *
     */
    public static RealVector getMVecForOULeaf (RealMatrix eMat, RealVector traitsValues, RealVector dVec) {
        return eMat.preMultiply(traitsValues).add(dVec);
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

    // (3) Lmr for internal node
    /*
     * For both BM and OU
     * L =
     * C - 0.25 * E * (A + L).inverse * E.transpose
     */
    public static RealMatrix getLMatForOUIntNode (RealMatrix cMat, RealMatrix eMat, RealMatrix eAPlusLInv) {
        return cMat.subtract(eAPlusLInv.multiply(eMat.transpose().scalarMultiply(0.25)));
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

    public static RealMatrix getPhiRM(Node aNode, RealMatrix phiMat) {

        double nodeBranchLength = aNode.getLength();

        return exp(phiMat.scalarMultiply(-nodeBranchLength));
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
    public static RealMatrix getOUVarianceRM (Node node, RealMatrix sigmaRM, RealMatrix sigmaeRM, RealMatrix pMat, RealMatrix inverseP, EigenDecomposition decompositionH, int nTraits) {
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

    public static RealMatrix getInverseVarianceRMForOU (Node node, double[] vcvMatDetArr, RealMatrix sigmaRM, RealMatrix sigmaeRM, RealMatrix pMat, RealMatrix inverseP, EigenDecomposition decompositionH, int nTraits) {

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
    public static void pruneOUPCM (Node node, int nTraits, List<RealVector> traitsValuesList,
                                   List<RealMatrix> lMatList, List<RealVector> mVecList, double[] rArr,
                                   RealMatrix sigmaRM, RealMatrix sigmaeRM, RealVector thetaVec, RealMatrix alphaMat,
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
            RealMatrix invVCVMat = getInverseVarianceRMForOU(child, vcvMatDetArr, sigmaRM, sigmaeRM, pMat, inverseP, decompositionH, nTraits);

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
                RealVector traitsVec = traitsValuesList.get(childIdx);

                // set the L matrix
                thisNodeLMat = thisNodeLMat.add(getLMatForOULeaf(cMat));

                // set r value
                double r = getRForOULeaf(aMat, traitsVec, bVec, f);
                thisNodeR += getRForOULeaf(aMat, traitsVec, bVec, f);

                // set m vector
                thisNodeMVec = thisNodeMVec.add(getMVecForOULeaf(eMat, traitsVec, dVec));
            } else {

                pruneOUPCM(child, nTraits, traitsValuesList, lMatList, mVecList, rArr, sigmaRM, sigmaeRM, thetaVec, alphaMat, pMat, inverseP, decompositionH, identity, vcvMatDetArr, negativeTwoAplusLDetArr);

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
}
