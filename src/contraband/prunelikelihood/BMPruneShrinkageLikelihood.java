package contraband.prunelikelihood;

import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.utils.PruneLikelihoodUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import outercore.parameter.KeyRealParameter;
import java.util.List;
import java.util.Random;

public class BMPruneShrinkageLikelihood extends PruneLikelihoodProcess {
    final public Input<Double> deltaInput = new Input<>("delta", "Shrinkage parameter for correlations, either sampled or given.");
    final public Input<Boolean> includePopVarInput = new Input<>("includePopVar", "if including population variance or not.", false);
    final public Input<RealParameter> popVarInput = new Input<>("popVar", "population variance.");
    final public Input<Double> deltaVarInput = new Input<>("deltaVar", "Shrinkage parameter for population variance, either sampled or given.");
    final public Input<KeyRealParameter> populationTraitsInput = new Input<>("populationTraits","Trait values for calculating the population noise.");

    private RealMatrix traitRM;
    private double delta;
    private double lambda;
    private double popVar;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // get delta: shrinkage parameter
        if (deltaInput.get() == null) {
            throw new RuntimeException("BMPruneShrinkageLikelihood::NodeMath is required for pmc likelihood.");
        }
        delta = deltaInput.get();

        // the real matrix that has the trait values for species
        traitRM = new Array2DRowRealMatrix(new double[getNSpecies()][getNTraits()]);
        PruneLikelihoodUtils.populateTraitValuesMatrix(traitsValuesInput.get(), treeInput.get(), getNTraits(), traitRM);

        if (includePopVarInput.get()) {
            /*
             * Consider population variance
             */
            setPopSE(true);
            if (populationTraitsInput.get()!=null) {
                lambda = deltaVarInput.get();
                // Estimate population variance from a sample
                RealMatrix popTraitsRM = populationTraitMatrix(populationTraitsInput.get());
                traitRM = PruneLikelihoodUtils.populateTraitValueMatrixEstimatedPopulationVariance(popTraitsRM, traitRM, getNTraits(), lambda);
                getNodeMath().estimateCorrelations(popTraitsRM);
            } else{
                getNodeMath().estimateCorrelations(traitRM);
                if (popVarInput.get() != null){
                    // Estimate population variance from input
                    popVar = popVarInput.get().getValue();
                    traitRM = PruneLikelihoodUtils.populateTraitValueMatrixGivenPopulationVariance(traitRM, popVar, getNTraits());
                } else {
                    lambda = deltaVarInput.get();
                    // Estimate population variance from data
                    traitRM = PruneLikelihoodUtils.populateTraitValueMatrixEstimatedPopulationVariance(traitRM, traitRM, getNTraits(), lambda);
                }
            }
        } else {
            /*
             * Ignore population variance
             */
            setPopSE(false);
            getNodeMath().estimateCorrelations(traitRM);
        }

        // linear shrinkage estimation of correlations
        // estimated in the beginning and fixed in the mcmc chain
        getNodeMath().populateShrinkageEstimation(delta);
        getNodeMath().populateTraitRateMatrix();
        getNodeMath().populateInverseTraitRateMatrix();
        getNodeMath().populateTransformedTraitValues(traitRM);
        setTraitValuesArr(getNodeMath().getTransformedTraitValues());
    }

    @Override
    public double calculateLogP() {
        // if the trait rate matrix has been changed
        // the data has to be transformed.
        if (getNodeMath().updateParameters()) {
            getNodeMath().populateTraitRateMatrix();
            getNodeMath().populateInverseTraitRateMatrix();
            getNodeMath().populateTransformedTraitValues(traitRM);
        }

        super.populateLogP();

        return getLogP();
    }

    private RealMatrix populationTraitMatrix(KeyRealParameter traitsValues){
        int nTraits = traitsValues.getMinorDimension1();
        int nSpecies = traitsValues.getMinorDimension2();
        RealMatrix traitRM = new Array2DRowRealMatrix(new double[nSpecies][nTraits]);
        String[] keys = traitsValues.getKeys();
        for (int i = 0; i < nSpecies; i ++) {
            // get all traits values for this species
            Double[] traitForSpecies = traitsValues.getRowValues(keys[i]);
            for (int j= 0; j < nTraits; j ++) {
                traitRM.setEntry(i, j, traitForSpecies[j]);
            }
        }
        return traitRM;
    }

    @Override
    protected void calculateLmrForTips(NodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        PruneLikelihoodUtils.populateLmrForTipWithShrinkage(nodeMath, traitValuesArr, nTraits, nodeIdx);
    }

    @Override
    protected void calculateLmrForInternalNodes(NodeMath nodeMath, int nTraits, int nodeIdx) {
        PruneLikelihoodUtils.populateLmrForInternalNodeWithShrinkage(nodeMath, nTraits, nodeIdx);
    }

    @Override
    protected double calculateLikelihood(NodeMath nodeMath, double l0, double[] m0, double r0, int rootIdx){
        double vCD = nodeMath.getVarianceForNode(rootIdx);
        double root2Subtract = -0.5 * nodeMath.getTraitRateMatrixDeterminant() - ((getNTraits() / 2.0) * Math.log(2 * Math.PI * vCD));
        return MatrixUtilsContra.vecTransScalarMultiply(nodeMath.getRootValuesArr(),
                l0, getNTraits()) +
                MatrixUtilsContra.vectorDotMultiply(nodeMath.getRootValuesArr(), m0) +
                r0 +
                nodeMath.getLikelihoodForSampledAncestors() -
                root2Subtract;

    }

    @Override
    public List<String> getArguments() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public List<String> getConditions() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        // TODO Auto-generated method stub

    }
}
