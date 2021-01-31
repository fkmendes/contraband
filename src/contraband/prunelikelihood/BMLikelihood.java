package contraband.prunelikelihood;

import contraband.math.GeneralNodeMath;
import contraband.math.MatrixUtilsContra;
import contraband.utils.MorphologyLikelihoodUtils;


public class BMLikelihood extends MorphologyLikelihood {


    private int nTraits;
    private boolean transformData;
    private boolean restrict;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        restrict = restrictInput.get();
        nTraits = traitInput.get().getTotalTraitNr();
        transformData = traitInput.get().getTransformDataFlag();
        nodeMathInput.get().updateSigmaMatrix();
        nodeMathInput.get().operateOnTraitRateMatrix();
        nodeMathInput.get().operateOnInvTraitRateMatrix();
        if(transformData) {
            traitInput.get().transformTraitData(nodeMathInput.get().getTraitRateMatrix());
        }

    }

    @Override
    protected boolean updateParameters() {
        // update BM model parameters
        boolean updateTraitRateMatrix = false;
        if(nodeMathInput.isDirty()){
            updateTraitRateMatrix  = nodeMathInput.get().updateSigmaMatrix();
        }

        return updateTraitRateMatrix;
    }

    @Override
    public double calculateLogP (){
        boolean update = updateParameters();

        if(update) {
            nodeMathInput.get().checkNearlySingularMatrix();

            if(nodeMathInput.get().getSingularMatrixFlag()){
                return Double.NEGATIVE_INFINITY;
            }

            nodeMathInput.get().operateOnTraitRateMatrix();
            nodeMathInput.get().operateOnInvTraitRateMatrix();
        }

        // update trait values (liabilities)
        traitInput.get().updateTraitValuesArr(update, nodeMathInput.get().getTraitRateMatrix());

        super.populateLogP();

        return getLogP();
    }


    @Override
    protected void populateAbCdEfForTips (GeneralNodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {
        if(nodeMath.getMatrixParamsFlag()){
            MorphologyLikelihoodUtils.populateACEfMatrixForTips(nodeMath, branchLength, nTraits, nodeIdx);
        } else {
            if(nodeMath.getPopVarianceFlag()) {
                double variance = branchLength * nodeMath.getTraitRate() + nodeMath.getPopVarianceMatrix()[0];
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
            if(nodeMath.getPopVarianceFlag()) {
                double variance = branchLength * nodeMath.getTraitRate();
                MorphologyLikelihoodUtils.populateACEf(nodeMath, variance, nTraits, nodeIdx);
            } else {
                MorphologyLikelihoodUtils.populateACEf(nodeMath, branchLength, nTraits, nodeIdx);
            }
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

        double root2Subtract = 0.0;
        if(restrict) {
            double vCD = nodeMath.getVarianceForNode(rootIdx);
            root2Subtract = -0.5 * nodeMath.getTraitRateMatrixDeterminant() - ((nTraits / 2.0) * Math.log(2 * Math.PI * vCD));
        }

        if(transformData){
            return MatrixUtilsContra.vecTransScalarMultiply(nodeMath.getRootValuesArr(),
                    l0[0], nTraits) +
                    MatrixUtilsContra.vectorDotMultiply(nodeMath.getRootValuesArr(), m0) +
                    r0 +
                    nodeMath.getLikelihoodForSampledAncestors() - root2Subtract;
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
                        r0 + nodeMath.getLikelihoodForSampledAncestors() - root2Subtract;
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
