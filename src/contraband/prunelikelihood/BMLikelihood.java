package contraband.prunelikelihood;

import contraband.math.GeneralNodeMath;
import contraband.math.MatrixUtilsContra;
import contraband.utils.MorphologyLikelihoodUtils;


public class BMLikelihood extends MorphologyLikelihood {


    private int nTraits;
    private boolean transformData;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        nTraits = traitInput.get().getTotalTraitNr();
        transformData = traitInput.get().getTransformDataFlag();
        nodeMathInput.get().updateSigmaMatrix();
        nodeMathInput.get().operateOnTraitRateMatrix();
        nodeMathInput.get().operateOnInvTraitRateMatrix();
    }

    @Override
    protected void updateParameters() {
        // update BM model parameters
        boolean updateSigmaMatrix = nodeMathInput.get().updateSigmaMatrix();
        if(updateSigmaMatrix) {
            nodeMathInput.get().operateOnTraitRateMatrix();
            nodeMathInput.get().operateOnInvTraitRateMatrix();
        }

        // update trait values (liabilities)
        traitInput.get().updateTraitValuesArr(updateSigmaMatrix);
    }

    @Override
    public double calculateLogP (){
        updateParameters();

        super.populateLogP();

        return getLogP();
    }


    @Override
    protected void populateAbCdEfForTips (GeneralNodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {
        if(nodeMath.getMatrixParamsFlag()){
            MorphologyLikelihoodUtils.populateACEfMatrixForTips(nodeMath, branchLength, nTraits, nodeIdx);
        } else {
            if(nodeMath.getPopVarianceFlag()) {
                double variance = branchLength + nodeMath.getPopVarianceMatrix()[0];
                MorphologyLikelihoodUtils.populateACEf(nodeMath, variance, nTraits, nodeIdx);
            } else {
                MorphologyLikelihoodUtils.populateACEf(nodeMath, branchLength, nTraits, nodeIdx);
            }
        }
    }

    @Override
    protected void populateAbCdEfForInternalNodes (GeneralNodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {
        if(nodeMath.getMatrixParamsFlag()){
            MorphologyLikelihoodUtils.populateACEfMatrixForIntNodes(nodeMath, branchLength, nTraits, nodeIdx);
        } else {
            MorphologyLikelihoodUtils.populateACEf(nodeMath, branchLength, nTraits, nodeIdx);
        }
    }

    @Override
    protected void populateLmrForTips(GeneralNodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        if(transformData) {
            MorphologyLikelihoodUtils.populateLmrForTipTransform(nodeMath, traitValuesArr, nTraits, nodeIdx);
        } else {
            if(nodeMath.getMatrixParamsFlag()) {
                MorphologyLikelihoodUtils.populateLmrMatrixForTip(nodeMath, traitValuesArr, nTraits, nodeIdx);
            } else {
                MorphologyLikelihoodUtils.populateLmrForTip(nodeMath, traitValuesArr, nTraits, nodeIdx);
            }
        }
    }

    @Override
    protected void populateLmrForInternalNodes(GeneralNodeMath nodeMath, int nTraits, int nodeIdx) {
        if(transformData) {
            MorphologyLikelihoodUtils.populateLmrForInternalNodeTransform(nodeMath, nTraits, nodeIdx);
        } else {
            if(nodeMath.getMatrixParamsFlag()) {
                MorphologyLikelihoodUtils.populateLmrMatrixForIntNode(nodeMath, nTraits, nodeIdx);
            } else {
                MorphologyLikelihoodUtils.populateLmrForIntNode(nodeMath, nTraits, nodeIdx);
            }
        }
    }

    @Override
    protected double calculateLikelihood(GeneralNodeMath nodeMath, double[] l0, double[] m0, double r0, int nTraits, int rootIdx){
        if(transformData){
            return MatrixUtilsContra.vecTransScalarMultiply(nodeMath.getRootValuesArr(),
                    l0[0], nTraits) +
                    MatrixUtilsContra.vectorDotMultiply(nodeMath.getRootValuesArr(), m0) +
                    r0 +
                    nodeMath.getLikelihoodForSampledAncestors();
        } else {
            if(nodeMath.getMatrixParamsFlag()){
                return MatrixUtilsContra.tVecDotMatrixDotVec(nodeMath.getRootValuesArr(), l0, nTraits)
                       + MatrixUtilsContra.vectorDotMultiply(
                       nodeMath.getRootValuesArr(), m0) +
                       r0 + nodeMath.getLikelihoodForSampledAncestors();
            } else {
                return l0[0] * MatrixUtilsContra.tVecDotMatrixDotVec(
                        nodeMath.getRootValuesArr(),
                        nodeMath.getTraitRateMatrixInverse(),
                        nTraits) +
                        MatrixUtilsContra.vectorDotMultiply(
                                nodeMath.getRootValuesArr(),
                                m0) +
                        r0 + nodeMath.getLikelihoodForSampledAncestors();
            }
        }
    }

    @Override
    protected double calculateLikelihoodForSA (GeneralNodeMath nodeMath, int childIdx) {
        if(transformData) {
            return MatrixUtilsContra.vecTransScalarMultiply(
                    nodeMath.getSampledAncestorTraitsVec(),
                    nodeMath.getLMatForNode(childIdx)[0],
                    nTraits) +
                    MatrixUtilsContra.vectorDotMultiply(
                            nodeMath.getSampledAncestorTraitsVec(),
                            nodeMath.getMVecForNode(childIdx)) +
                    nodeMath.getRForNode(childIdx);
        } else {
            return nodeMath.getLMatForNode(childIdx)[0] *
                    MatrixUtilsContra.tVecDotMatrixDotVec(
                            nodeMath.getSampledAncestorTraitsVec(),
                            nodeMath.getTraitRateMatrixInverse(),
                            nTraits) +
                    MatrixUtilsContra.vectorDotMultiply(
                            nodeMath.getSampledAncestorTraitsVec(),
                            nodeMath.getMVecForNode(childIdx)) +
                    nodeMath.getRForNode(childIdx);
        }
    }
}
