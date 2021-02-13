package testdrivers;

import contraband.math.MatrixUtilsContra;

public class OUModelParamsTestDriver2 {
    public static void main(String[] args) {

        // (1) initialize the test
        int nTraits = 4;
        double[] aPlusLInverse = new double[nTraits * nTraits];
        double logDetNegativeTwoAplusL = 1.0;

        double[] bPlusM = new double[nTraits];
        double[] mVec = new double[nTraits];
        double[] dToSubtract = new double[nTraits];

        //r = f + r + 0.5 * nTraits * log(2*PI) - 0.5 * log(det(-2 * (A + L))) - 0.25 * (b + m).transpose * (A + L).inverse * (b + m)
        //mVec = bVec.add(mVec);
        //r = r + f + 0.5 * nTraits* LOGTWOPI - 0.5 * logDetVNode  - 0.25 * mVec.dotProduct(invAPlusLRM.preMultiply(mVec));
        MatrixUtilsContra.vectorAdd(nodeMath.getBVecForNode(nodeIdx), nodeMath.getMVecForNode(nodeIdx), bPlusM);
        double r = nodeMath.getRForNode(nodeIdx) + nodeMath.getfForNode(nodeIdx)
                + 0.5 * nTraits* LOGTWOPI
                - 0.5 * nodeMath.getLogDetNegativeTwoAplusL()
                - 0.25 * MatrixUtilsContra.tVecDotMatrixDotVec(bPlusM, nodeMath.getAPlusLInverse(), nTraits);
        nodeMath.setRForNode(nodeIdx, r);

        //m = d - 0.5 * E * (A + L).inverse * （bVec + mVec)
        //mVec = bVec.add(mVec);
        //dVec.subtract(eAPlusLInv.transpose().preMultiply(mVec).mapMultiply(0.5));
        // m_non-tip = - 0.5 *  E * (A + L).inverse * m_i
        double[] eAPlusLInv = new double[nTraits * nTraits];
        double[] eAplusLInvTrans = new double[nTraits * nTraits];
        MatrixUtilsContra.matrixMultiply(nodeMath.getEMatForNode(nodeIdx), nodeMath.getAPlusLInverse(), nTraits, nTraits, eAplusLInvTrans);
        MatrixUtilsContra.matrixTranspose(eAplusLInvTrans, nTraits, eAPlusLInv);
        MatrixUtilsContra.matrixTransPreMapMultiply(bPlusM, eAPlusLInv, -0.5, nTraits, nTraits, dToSubtract);
        MatrixUtilsContra.vectorSubtract(nodeMath.getDVecForNode(nodeIdx), dToSubtract, mVec);
        nodeMath.setMVecForNode(nodeIdx, mVec);

        //C - 0.25 * E * (A + L).inverse * E.transpose
        //cMat.subtract(eAPlusLInv.multiply(eMat.transpose()).scalarMultiply(0.25));
        // L = C - 0.25 * E * (A + L).inverse * E.transpose
        double[] eMatTransScalar = new double [nTraits * nTraits];
        double[] cToSubtract = new double [nTraits * nTraits];
        double[] lMat = new double [nTraits * nTraits];

        MatrixUtilsContra.matrixTransScalar(nodeMath.getEMatForNode(nodeIdx), 0.25, nTraits, eMatTransScalar);
        MatrixUtilsContra.matrixMultiply(eAPlusLInv, eMatTransScalar, nTraits, nTraits, cToSubtract);
        MatrixUtilsContra.matrixSubtract(nodeMath.getCMatForNode(nodeIdx), cToSubtract, nTraits, lMat);
        nodeMath.setLForNode(nodeIdx, lMat);
    }
}
