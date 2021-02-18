package contraband.prunelikelihood;

import beast.core.Input;
import contraband.math.GeneralNodeMath;
import contraband.math.MatrixUtilsContra;
import contraband.utils.MorphologyLikelihoodUtils;

public class OULikelihood extends MorphologyLikelihood {
    final public Input<OUModelParameter> modelParameterInput = new Input<>("params","Parameters for OU model", Input.Validate.REQUIRED);

    OUModelParameter modelParameter;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        nodeMathInput.get().intOUModelParams();
        nodeMathInput.get().updateSigmaMatrix();
        modelParameter = modelParameterInput.get();
        modelParameter.updateOUModelParams();
        modelParameter.performAlphaEigenDecomposition();
    }

    @Override
    protected boolean updateParameters() {
        // update OU model parameters
        boolean updateTraitRateMatrix = false;
        if(nodeMathInput.isDirty()){
            updateTraitRateMatrix  = nodeMathInput.get().updateSigmaMatrix();
        }

        if(modelParameterInput.isDirty()) {
            modelParameter.updateOUModelParams();
            //modelParameter.performAlphaEigenDecomposition();
        }

        return updateTraitRateMatrix;
    }

    @Override
    public double calculateLogP() {
        boolean update = updateParameters();

        if(modelParameter.isSingularAlphaMatrix()) { return Double.NEGATIVE_INFINITY; }

        // update trait values (liabilities)
        traitInput.get().updateTraitValuesArr(update, nodeMathInput.get().getTraitRateMatrix());

        super.populateLogP();

        return getLogP();
    }

    @Override
    protected void populateAbCdEfForTips(GeneralNodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {
        MorphologyLikelihoodUtils.populatePhiMatrixForNode(traitInput.get().getTotalTraitNr(), branchLength, modelParameter.getAlphaMat(), nodeMath);
        MorphologyLikelihoodUtils.populateOmegaVecForNode(nodeIdx, traitInput.get().getTotalTraitNr(), nodeMath.getPhiMatForNode(nodeIdx), modelParameter, nodeMath);
        nodeMath.populateOUVarianceMatrixForNode(nodeIdx, branchLength, modelParameter, true);
        nodeMath.operateOnVarianceMatrix();
        MorphologyLikelihoodUtils.populateAbCdEf(nodeMath, nTraits, nodeIdx);
    }

    @Override
    protected void populateAbCdEfForInternalNodes (GeneralNodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {
        MorphologyLikelihoodUtils.populatePhiMatrixForNode(traitInput.get().getTotalTraitNr(), branchLength, modelParameter.getAlphaMat(), nodeMath);
        MorphologyLikelihoodUtils.populateOmegaVecForNode(nodeIdx, traitInput.get().getTotalTraitNr(), nodeMath.getPhiMatForNode(nodeIdx), modelParameter, nodeMath);
        nodeMath.populateOUVarianceMatrixForNode(nodeIdx, branchLength, modelParameter, false);
        nodeMath.operateOnVarianceMatrix();
        MorphologyLikelihoodUtils.populateAbCdEf(nodeMath, nTraits, nodeIdx);
    }

    @Override
    protected void populateLmrForTips(GeneralNodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        MorphologyLikelihoodUtils.populateLmrForOUTip(nodeMath, traitValuesArr, nTraits, nodeIdx);
    }

    @Override
    protected void populateLmrForInternalNodes(GeneralNodeMath nodeMath, int nTraits, int nodeIdx) {
        MorphologyLikelihoodUtils.populateLmrForOUIntNode(nodeMath, nTraits, nodeIdx);
    }

    @Override
    protected double calculateLikelihood(GeneralNodeMath nodeMath, double[] l0, double[] m0, double r0, int nTraits, int rootIdx) {
        return MatrixUtilsContra.tVecDotMatrixDotVec(nodeMath.getRootValuesArr(), l0, nTraits)
                + MatrixUtilsContra.vectorDotMultiply(
                nodeMath.getRootValuesArr(), m0) +
                r0 + nodeMath.getLikelihoodForSampledAncestors();
    }

    @Override
    protected double calculateLikelihoodForSA (GeneralNodeMath nodeMath, int childIdx) {
        return MatrixUtilsContra.tVecDotMatrixDotVec(
                        nodeMath.getSampledAncestorTraitsVec(),
                        nodeMath.getLMatForNode(childIdx),
                        traitInput.get().getTotalTraitNr()) +
                MatrixUtilsContra.vectorDotMultiply(
                        nodeMath.getSampledAncestorTraitsVec(),
                        nodeMath.getMVecForNode(childIdx)) +
                nodeMath.getRForNode(childIdx);
    }

}

