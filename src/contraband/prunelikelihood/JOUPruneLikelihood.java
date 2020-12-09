package contraband.prunelikelihood;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.*;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import contraband.valuewrappers.OneValueContTraits;
import org.apache.commons.math3.linear.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

@Description("This class implements likelihood for continuous traits under Ornstein–Uhlenbeck process.\n" +
        "The calculation uses Venelin's PCM likelihood.")


public class JOUPruneLikelihood extends OUPruneLikelihoodProcess {


    @Override
    public void initAndValidate() {
        super.initAndValidate();
        nodeMath.initValidateJOU();
    }

    @Override
    public double calculateLogP() {
        nodeMath.updateJOUParameters();
        nodeMath.populateAlphaMatrix();
        nodeMath.performAlphaDecomposition(nodeMath.getAlphaMatrix());

        super.populateLogP();

        return getLogP();
    }

    @Override
    protected void populateAbCdEf (OUNodeMath nodeMath, Node node, int nodeIdx) {
        // For OU, variance matrix, Phi and Omega need to be calculated for this node.
        RealMatrix phiRM = calculatePhiMatrix(node, nodeMath);
        RealVector omegaVec = calculateOmegaVector(node, nodeMath, phiRM);
        nodeMath.setPhiMatForNode(nodeIdx, phiRM);

        // inverse of variance-covariance of this node
        nodeMath.populateVarianceCovarianceMatrix(node);
        nodeMath.populateVarianceCovarianceMatrixForJOU(nodeIdx);
        nodeMath.setVarianceCovarianceMatrix(node, nodeIdx);

        RealMatrix aMat = OUPruneUtils.getAMatForOU(nodeMath.getInverseVarianceMatrix());
        RealMatrix eMat = OUPruneUtils.getEMatForOU(phiRM, nodeMath.getInverseVarianceMatrix());
        RealMatrix cMat = OUPruneUtils.getCMatForOU(phiRM, eMat);
        double f = OUPruneUtils.getFforOU(omegaVec, nodeMath.getInverseVarianceMatrix(), nodeMath.getVCVMatDet(), phiRM.getColumnDimension());
        RealVector bVec = OUPruneUtils.getBVecForOU(nodeMath.getInverseVarianceMatrix(), omegaVec);
        RealVector dVec = OUPruneUtils.getDVecForOU(eMat, omegaVec);

        nodeMath.setAbCdEfOmegaPhiForNode(nodeIdx, aMat, bVec, cMat, dVec, eMat, f, omegaVec, phiRM);
    }

    /*
     * The method calculates omega vector for OU model with jumps
     * H <- alpha matrix
     * t <- branch length in time
     * mj <- jump mean vector
     * phiRM <- expm(-t*H)
     * I <- diag(nTraits)
     * if (a jump happens to a node) {
     *  omega <- phiRM%*%mj +  (I-phiRM)%*%Theta
     * } else {
     *  omega <- (I-phiRM)%*%Theta
     * }
     */
    @Override
    protected RealVector calculateOmegaVector (Node node, OUNodeMath nodeMath, RealMatrix phiRM) {
        int xi = nodeMath.getJumpForNode(node.getNr());

        RealVector iMinusPhiTheta = OUPruneUtils.getOmegaVec(nodeMath.getThetaForNode(node.getNr()), phiRM, nodeMath.getIdentityMatrix());

        if(xi == 0) {
            return iMinusPhiTheta;
        } else{
            return  OUPruneUtils.getOmegaVec(nodeMath.getJumpMeanVec(), phiRM, iMinusPhiTheta);
        }
    }

}
