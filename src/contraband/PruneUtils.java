package contraband;

import beast.evolution.tree.Node;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import java.util.List;

public class PruneUtils {

    static double LOGTWOMATHPI = Math.log(2 * Math.PI);

    /*
     * This sets the evolutionary rate matrix in place for a node, under BM.
     *
     * This is a matrix instead of scalar because we can have
     * from 1 up to 'k' traits.
     *
     * Under BM, note that we do not have to worry about the inverse of eigen
     * matrix 'P', invP, and its transpose (all of these consist of the identity matrix),
     * so we can skip the operation invP * vcvMatForBranch * invP^T.
     */
    public static void setVCVMatForBranchBM(Node aNode, RealMatrix evolRateMat, List<RealMatrix> vcvMatList) {
        double t = aNode.getLength(); // length of the branch subtending this node
        vcvMatList.set(aNode.getNr(), evolRateMat.multiply(evolRateMat.transpose()).scalarMultiply(t)); // make evol rate matrix symmetric
    }

    /*
     * Matrix 'A' is -0.5 of the inverse of vcvMatForBranch.
     *
     * There will be one aMat per node.
     */
    public static void setAMat(Node aNode, RealMatrix invVCVMatForBranch, List<RealMatrix> aMatList) {
        aMatList.set(aNode.getNr(), invVCVMatForBranch.scalarMultiply(-0.5));
    }

    /*
     * The 'E' matrix in BM is equivalent to invVCVMatForBranch,
     * because matrix 'tPhiMat' is the identity matrix.
     *
     * Under OU, then 'tPhiMat' is the tranpose of phi, which in turn represents ...
     * and is/is not a parameter that is used for ... [Rong: FILL THIS OUT]
     *
     * There will be one eMat per node.
     */
    public static void setEMatOU(Node aNode, RealMatrix invVCVMatForBranch, RealMatrix tPhiMat, List<RealMatrix> eMatList) {
        eMatList.set(aNode.getNr(), tPhiMat.multiply(invVCVMatForBranch));
    }

    /*
     * The 'B' vector (nTrait size) in BM is a vector of 0's
     * because the omega vector 'omegaVec' is also a vector of 0's.
     *
     * Under OU, the 'omegaVec' vector represents ... and is/is not a parameter
     * that is used for ... [Rong: FILL THIS OUT]
     *
     * There will be one bVec per node.
     */
    public static void getBVecOU(RealMatrix invVCVMatForBranch, RealVector omegaVec, RealVector bVec) {
        ;
    }

    /*
     * The 'C' matrix, 'cMat', is equivalent to the 'A' matrix under BM.
     *
     * There will be one cMat per node.
     */
    public static void setCMatOU(Node aNode, RealMatrix eMat, RealMatrix phiMat, List<RealMatrix> cMatList) {
        cMatList.set(aNode.getNr(), eMat.multiply(phiMat).scalarMultiply(-0.5));
    }

    /*
     * Under BM, the 'D' vector, 'dVec', will be the product of
     * -invVCVMatForBranch * omegaVec (which under BM is a vector
     * of zeros). So dVec for BM is a vector of 0's.
     *
     * There will be one eMat per node.
     */
    public static void setDVecOU(Node aNode, RealMatrix eMat, RealVector omegaVec, List<RealVector> dVecList) {
        dVecList.set(aNode.getNr(), eMat.scalarMultiply(-1.0).preMultiply(omegaVec));
    }

    /*
     * Under BM, omegaVec is a vector of 0's, and so we can
     * ignore the first term in the equation that defines 'f'
     * in Mitov et al. (2019).
     *
     * 'f' is then just a function of the log-determinant of
     * the vcvMatForBranch.
     *
     * There will be one f per node.
     */
    public static void setF(Node aNode, int nTraits, double vcvMatDet, double[] fArr) {
        fArr[aNode.getNr()] = -0.5 * nTraits * LOGTWOMATHPI + Math.log(vcvMatDet);
    }

    public static RealMatrix getLMatForLeaf(RealMatrix cMat) {
        RealMatrix tCMat = cMat.transpose();
        RealMatrix lMat = cMat.add(tCMat).scalarMultiply(0.5);

        return lMat;
    }

    /*
     * Document here
     */
    public static void setLMatForIntNode(Node aNode, List<RealMatrix> lMatList) {
        int thisNodeIdx = aNode.getNr();
        RealMatrix thisNodeLMat = lMatList.get(thisNodeIdx);

        List<Node> children = aNode.getChildren();

        // running sum
        for (Node child : children) {
            int childIdx = child.getNr();
            thisNodeLMat = thisNodeLMat.add(lMatList.get(childIdx));
        }

        lMatList.set(thisNodeIdx, thisNodeLMat);
    }

    /*
     * Under BM, bVec is a vector of 0's, so we can ignore the second
     * term inside the summation in the equation that defines r in
     * Mitov et al. (2019).
     *
     * There will be one r per leaf (and one r per internal node).
     */
    public static double getRForLeafBM(RealMatrix aMat, RealVector traitsValues, RealVector tTraitsValues, double f) {
        double r = aMat.preMultiply(traitsValues).dotProduct(traitsValues);

        return r;
    }

    /*
     * Under BM, dVec is a vector of 0's, se we can ignore the first
     * term inside the summation in the equation that defines mBec in
     * Mitov et al. (2019).
     *
     * There will be one eMat per leaf (and one eMat per internal node).
     */
    public static RealVector getMVecForLeafBM(RealMatrix eMat, RealVector traitsValues) {
        RealVector mVec = eMat.preMultiply(traitsValues);

        return mVec;
    }


}
