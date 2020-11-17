package contraband.prunelikelihood;

import beast.evolution.tree.Node;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

public class OUPruneUtils {
    private static double LOGTWOPI = FastMath.log(2 * Math.PI);

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

    // (3) trait rate matrix and branch length in time
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
                    f = 1 - Math.exp(-lambda * nodeBranchLength) / lambda;
                    //fLambdaP_1SigmaP_t[i][j] = f * P_1SigmaP_t.getEntry(i,j);
                    variance.setEntry(i, j, f * variance.getEntry(i, j));
                }
            }
        }

        // variance matrix = P * (fLambda * P_1SigmaP_t) * t(P)
        variance = pMat.multiply(variance).multiply(pMat.transpose());


        if (node.isLeaf() && sigmaeRM!=null) {
            variance = variance.add(sigmaRM.multiply(sigmaRM.transpose()));
        }

        return variance;
    }
}
