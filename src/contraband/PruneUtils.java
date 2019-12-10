package contraband;

import beast.evolution.tree.Node;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class PruneUtils {

    /*
     * This returns the evolutionary rate matrix in place under BM
     *
     * This is a matrix instead of scalar because we can have
     * from 1 up to 'k' traits
     *
     * Under BM, note that we do not have to worry about the inverse of eigen
     * matrix 'P', invP, and its transpose (all of these consist of the identity matrix),
     * so we can skip the operation invP * vcvMatForBranch * invP^T
     */
    public static RealMatrix getVCVMatForBranchInPlaceBM(Node aNode, RealMatrix evolRateMat) {
        double t = aNode.getLength(); // length of the branch subtending this node
        RealMatrix vcvMatForBranch = evolRateMat.multiply(evolRateMat.transpose()); // make evol rate matrix symmetric
        return vcvMatForBranch.scalarMultiply(t);
    }

    /*
     * Matrix 'A' is -0.5 of the inverse of vcvMatForBranch
     */
    public static RealMatrix getAMatInPlace(RealMatrix invVCVMatForBranch) {
        RealMatrix aMat = invVCVMatForBranch.scalarMultiply(-0.5);

        return aMat;
    }

    /*
     * The 'E' matrix in BM is equivalent to invVCVMatForBranch,
     * because matrix 'tPhiMat' is the identity matrix.
     *
     * Under OU, then 'tPhiMat' is the tranpose of phi, which in turn represents ...
     * and is/is not a parameter that is used for ... [Rong: FILL THIS OUT]
     */
    public static RealMatrix getEMatOU(RealMatrix invVCVMatForBranch, RealMatrix tPhiMat) {
        RealMatrix eMat = tPhiMat.multiply(invVCVMatForBranch);

        return eMat;
    }

    /*
     * The 'B' vector (nTrait size) in BM is a vector of 0's
     * because the omega vector 'omegaVec' is also a vector of 0's
     *
     * Under OU, the 'omegaVec' vector represents ... and is/is not a parameter
     * that is used for ... [Rong: FILL THIS OUT]
     */
    public static void getBVecOU(RealMatrix invVCVMatForBranch, RealVector omegaVec, RealVector bVec) {
        ;
    }

    /*
     * The 'C' matrix, 'cMat', is equivalent to the 'A' matrix under BM.
     */
    public static RealMatrix getCMatOU(RealMatrix eMat, RealMatrix phiMat) {
        RealMatrix cMat = eMat.multiply(phiMat).scalarMultiply(-0.5);

        return cMat;
    }

    /*
     * Under BM, the 'D' vector, 'dVec', will be the product of
     * -invVCVMatForBranch * omegaVec (which under BM is a vector
     * of zeros). So dVec for BM is a vector of 0's
     */
    public static RealVector getDVecOU(RealMatrix eMat, RealVector omegaVec) {
        return eMat.scalarMultiply(-1.0).preMultiply(omegaVec);
    }

    /*
     * Under BM, f will be the product of
     * -0.5 and sumDetLog2PI (which under BM sumDetLog2PI is
     * nTraits * Math.log(2 * Math.PI) and Math.log(vCVMatDet)).
     *
     */
    public static double getf(int nTraits, double vCVMatDet) {
        double sumDetLog2PI = nTraits * Math.log(2 * Math.PI) + Math.log(vCVMatDet);
        return -0.5 * sumDetLog2PI;
    }

}
