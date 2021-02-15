package contraband.prunelikelihood;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import contraband.math.LUDecompositionForArray;
import contraband.math.MatrixUtilsContra;
import contraband.utils.MorphologyLikelihoodUtils;
import contraband.utils.NodeMathUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import outercore.parameter.KeyRealParameter;


// This class deals with SigmaX matrix in PCM likelihood.
// It can be used to populate either Sigma (trait rate matrix, i.e. evolutionary rate in unit variance),
// or Sigmae (variance covariance matrix for intraspecific variation)
public class SigmaMatrix extends CalculationNode {
    final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Evolutionary rates of traits. Diagonal elements in rate matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> rhoInput = new Input<>("correlation", "Correlations of traits. Off-diagonal elements in rate matrix.");
    final public Input<MorphologicalData> traitInput = new Input<>("trait","Morphological data set.");
    final public Input<Boolean> oneRateOnlyInput = new Input<>("oneRateOnly", "TRUE, if all traits share one evolutionary rate.", false);
    final public Input<Boolean> useUpperMatrixInput = new Input<>("upperMatrix", "TRUE, if sigmasq and correlations are from upper matrix", false);
    final public Input<RealParameter> covarianceInput = new Input<>("covariance", "cov_ij = sigma_i * sigma_j * rho_ij.");
    final public Input<RealParameter> inverseRhoInput = new Input<>("inverseMatrix", "Inverse of trait correlation matrix when using Gibbs sampling.");
    final public Input<Boolean> useShrinkageInput = new Input<>("shrinkage", "TRUE, if shrinkage method is used to estimate the trait correlations.", false);
    final public Input<KeyRealParameter> populationInput = new Input<>("population","Trait values for calculating the population noise. Continuous data only.");
    final public Input<RealParameter> deltaInput = new Input<>("delta", "Shrinkage intensity parameter.");

    private double sigmaValue;
    private double[] sigmaValues;
    private double[] rhoValues;
    private double[] inverseRhoMatrix;
    private double[] covarianceValues;

    private boolean useUpperMatrix;
    private boolean oneRateOnly;
    private boolean coEstimate;
    private boolean inverseMatrix;
    private boolean useShrinkage;

    private int nTraits;
    private int matDim;

    private double[] sigmaMatrix;
    private double[] storedSigmaMatrix;
    private double[] rhoMatrix;

    private double [] upperMatrix;
    private double [] transUpperMatrix;

    // using shrinkage
    private double [] unbiasedRho;
    private double delta;
    private RealMatrix populationTraitMatrix;

    // store and restore
    private double storedSigmaValue;
    private double[] storedSigmaValues;
    private double[] storedRhoValues;
    private double[] storedInverseRhoMatrix;
    private double[] storedRhoMatrix;
    private double[] storedCovariance;

    private double[] lu;
    private int[] pivot;
    private boolean[] evenSingular;
    private double [] identityMatrix;
    private double detInvRhoMatrix;
    private boolean singularMatrix;

    @Override
    public void initAndValidate() {
        useUpperMatrix = useUpperMatrixInput.get();
        oneRateOnly = oneRateOnlyInput.get();
        coEstimate = covarianceInput.get() != null;
        inverseMatrix = inverseRhoInput.get() != null;
        useShrinkage = useShrinkageInput.get();

        nTraits = traitInput.get().getTotalTraitNr();
        initArrays();
        initUDecomposition();

        sigmaMatrix = new double[nTraits * nTraits];
        storedSigmaMatrix = new double[nTraits * nTraits];
        rhoMatrix = new double[nTraits * nTraits];
        for (int i = 0; i < nTraits; i++){
            rhoMatrix[i * nTraits + i] = 1.0;
        }
        storedRhoMatrix = new double[nTraits * nTraits];

        // get elements on the diagonal
        if(oneRateOnly) {
            sigmaValue = sigmasqInput.get().getValue();
            if(sigmasqInput.get().getDimension() != 1) {
                sigmasqInput.get().setDimension(1);
            }
        } else {
            sigmaValues = new double[nTraits];
            storedSigmaValues = new double[nTraits];
            if(sigmasqInput.get().getDimension() != nTraits) {
                sigmasqInput.get().setDimension(nTraits);
            }
            sigmaValues = sigmasqInput.get().getDoubleValues();
        }

        // get elements on the off-diagonal
        if(inverseMatrix) {
            // if inverse, it is a full matrix
            rhoValues = new double[nTraits * nTraits];
            inverseRhoMatrix = new double[nTraits * nTraits];
            storedRhoValues = new double[nTraits * nTraits];
            storedInverseRhoMatrix = new double[nTraits * nTraits];
            initInverseRho();
            populateRhoValuesFromInverseMatrix();
        } else {
            // if not inverse, it is an upper diagonal matrix
            rhoValues = new double[nTraits * (nTraits - 1) / 2];
            storedRhoValues = new double[nTraits * (nTraits - 1) / 2];

            // use shrinkage method to estimate rho matrix from a population sample of a species
            if(useShrinkage) {
                populationTraitMatrix = MorphologyLikelihoodUtils.populateTraitMatrixForPopulationSample(populationInput.get());

                unbiasedRho = new double[nTraits * (nTraits - 1) / 2];

                estimateCorrelationUsingShrinkage();
            }

            else if(coEstimate) {
                // an upper diagonal matrix
                covarianceValues = new double[nTraits * (nTraits - 1) / 2];
                storedCovariance = new double[nTraits * (nTraits - 1) / 2];
                if (covarianceInput.get().getDimension() != (nTraits * (nTraits - 1) / 2)) {
                    covarianceInput.get().setDimension((nTraits * (nTraits - 1) / 2));
                }
                covarianceValues = covarianceInput.get().getDoubleValues();
            }
            
            else if (rhoInput.get()!=null){
                if (rhoInput.get().getDimension() != (nTraits * (nTraits - 1) / 2)) {
                    rhoInput.get().setDimension((nTraits * (nTraits - 1) / 2));
                }
                rhoValues = rhoInput.get().getDoubleValues();
            }
        }
    }

    private void initArrays() {
        matDim = nTraits * nTraits;
        upperMatrix = new double[matDim];
        transUpperMatrix = new double[matDim];
    }

    private void initUDecomposition() {
        lu = new double [matDim];
        pivot = new int[nTraits];
        evenSingular = new boolean[2];
        identityMatrix = new double[matDim];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, nTraits);
        }
    }

    private void initInverseRho(){
        if(inverseRhoInput.get().getDimension() != matDim) {
            inverseRhoInput.get().setDimension(matDim);
            // initiate an identity matrix
            for (int i = 0; i < matDim; i++){
                if(i % (nTraits + 1) == 0) {
                    inverseRhoInput.get().setValue(i, 1.0);
                    inverseRhoMatrix[i] = 1.0;
                } else {
                    inverseRhoInput.get().setValue(i, 0.0);
                    inverseRhoMatrix[i] = 0.0;
                }
            }
        } else {
            // get from input
            inverseRhoMatrix = inverseRhoInput.get().getDoubleValues();
        }
    }

    // getters
    public double [] getSigmaMatrix () { return sigmaMatrix; }

    public double [] getRhoMatrix () { return rhoMatrix; }

    public boolean getOneRateOnly () { return oneRateOnly; }

    // main
    public void populateSigmaValue () {
        sigmaMatrix[0] = sigmaValue;
    }

    public void populateRhoMatrix () {
        int m = 0;
        for (int i = 0; i < nTraits; i++) {
            for (int j = (i + 1); j < nTraits; j++){
                double cor = rhoValues[m];
                rhoMatrix[i * nTraits + j] = cor;
                rhoMatrix[j * nTraits + i] = cor;
                m++;
            }
        }
    }

    public void populateSigmaMatrix() {
        if(oneRateOnly) {
            if(coEstimate){
                populateSigmaForOneRateCo();
            } else {
                populateSigmaForOneRate();
            }
        } else {
            if(coEstimate){
                populateSigmaForMultipleRatesCo ();
            } else {
                populateSigmaForMultipleRates();
            }
        }
    }

    // sigma_ij = sigmasq * rho_ij
    private void populateSigmaForOneRate () {
        if(useUpperMatrix) {
            // sigma = t(upper) %*% upper
            NodeMathUtils.populateTraitRateMatrix(sigmaValue, rhoValues, upperMatrix, transUpperMatrix, nTraits, sigmaMatrix);
        } else {
            if(inverseMatrix) {
                for(int i = 0; i < matDim; i ++){
                    sigmaMatrix[i] = sigmaValue * rhoValues[i];
                }
            } else {
                NodeMathUtils.populateTraitRateMatrixDirectly(sigmaValue, rhoValues, nTraits, sigmaMatrix);
            }
        }
    }

    // sigma_i = sigmasq; sigma_ij = cov_ij
    private void populateSigmaForOneRateCo () {
        if(useUpperMatrix) {
            populateSigmaMatrixCo(sigmaValue, covarianceValues, upperMatrix, transUpperMatrix, nTraits, sigmaMatrix);
        } else {
            populateSigmaMatrixDirectlyCo(sigmaValue, covarianceValues, nTraits, sigmaMatrix);
        }
    }

    // sigma_ij =  sqrt(sigmasq_i) * sqrt(sigmasq_j) * rho_ij
    private void populateSigmaForMultipleRates () {
        if (useUpperMatrix) {
            // sigma = t(upper) %*% upper
            NodeMathUtils.populateTraitRateMatrix(sigmaValues, rhoValues, upperMatrix, transUpperMatrix, nTraits, sigmaMatrix, coEstimate);
        } else {
            if(inverseMatrix) {
                for(int i = 0; i < nTraits; i++){
                    for(int j = 0; j < nTraits; j++){
                        sigmaMatrix[i * nTraits + j] = FastMath.sqrt(sigmaValues[i]) * FastMath.sqrt(sigmaValues[j]) * rhoValues[i * nTraits + j];
                    }
                }
            } else {
                NodeMathUtils.populateTraitRateMatrixDirectly(sigmaValues, rhoValues, nTraits, sigmaMatrix);
            }
        }
    }

    // sigma_i = sigmasq_i; sigma_ij = cov_ij
    private void populateSigmaForMultipleRatesCo () {
        if (useUpperMatrix) {
            populateSigmaMatrixCo(sigmaValues, covarianceValues, upperMatrix, transUpperMatrix, nTraits, sigmaMatrix);
        } else {
            populateSigmaMatrixDirectlyCo(sigmaValues, covarianceValues, nTraits, sigmaMatrix);
        }
    }

    public boolean updateParameters(boolean sigmaValueType){
        boolean updateRho = false;
        boolean updateSigma = false;

        // update trait correlations
        if(rhoInput.isDirty()) {
            rhoValues = rhoInput.get().getDoubleValues();
            updateRho = true;
        }
        if(inverseRhoInput.isDirty()) {
            inverseRhoMatrix = inverseRhoInput.get().getDoubleValues();
            updateRho = true;
        }
        // update trait correlations
        if(covarianceInput.isDirty()) {
            covarianceValues = covarianceInput.get().getDoubleValues();
            updateRho = true;
        }

        // update trait evolutionary rates
        if(sigmasqInput.isDirty()) {
            if (oneRateOnly) {
                sigmaValue = sigmasqInput.get().getValue();
            } else {
                sigmaValues = sigmasqInput.get().getDoubleValues();
            }
            updateSigma = true;
       }

        boolean updateSigmaMatrix = updateRho || updateSigma;

        if(updateSigmaMatrix){
                if(inverseMatrix) {
                    populateRhoValuesFromInverseMatrix();
                }
                //TO DO: confirm the condition of updating matrix!!!
                //if(!sigmaValueType) {
                if(!oneRateOnly) {
                    // update the sigma matrix using rho values and sigma values
                    populateSigmaMatrix();
                    }
                else {
                    // if shared rho between trait rate matrix and population variance
                    // we only update the shared rate/variance among traits and the rho matrix
                    populateSigmaValue();
                    if (rhoInput.get() != null) {
                        populateRhoMatrix();
                    }
                }
        }
        return updateSigmaMatrix;
    }

    private void estimateCorrelationUsingShrinkage() {
        MorphologyLikelihoodUtils.populateUnbiasedRho(populationTraitMatrix, unbiasedRho);

        // R = delta * I + (1 - delta) * R^
        delta = deltaInput.get().getValue();
        MatrixUtilsContra.vectorMapMultiply(unbiasedRho, 1 - delta, rhoValues);
    }

    public void populateRhoValuesFromInverseMatrix() {
        // to get the original rho matrix because it will be used when calculating r parameter
        LUDecompositionForArray.ArrayLUDecomposition(inverseRhoMatrix, lu, pivot, evenSingular, nTraits);
        detInvRhoMatrix = LUDecompositionForArray.getDeterminant(lu, nTraits, evenSingular);
        try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evenSingular[1], nTraits, rhoValues);
        } catch (RuntimeException e) {
            singularMatrix = true;
        }
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }
    @Override
    public void store() {
        super.store();

        if(oneRateOnly){
            storedSigmaValue = sigmaValue;
        } else {
            System.arraycopy(sigmaValues, 0, storedSigmaValues, 0, nTraits);
        }

        System.arraycopy(rhoValues, 0, storedRhoValues, 0, rhoValues.length);
        if(inverseMatrix) {
            System.arraycopy(inverseRhoMatrix, 0, storedInverseRhoMatrix, 0, matDim);
        }

        System.arraycopy(sigmaMatrix, 0, storedSigmaMatrix, 0, sigmaMatrix.length);
        System.arraycopy(rhoMatrix, 0, storedRhoMatrix, 0, storedRhoMatrix.length);

        if(coEstimate){
            System.arraycopy(covarianceValues, 0, storedCovariance, 0, covarianceValues.length);
        }
    }

    @Override
    public void restore() {
        super.restore();
        if(oneRateOnly) {
            double tempSigmaValue = sigmaValue;
            sigmaValue = storedSigmaValue;
            storedSigmaValue = tempSigmaValue;
        } else {
            double[] tempSigmaValues = sigmaValues;
            sigmaValues = storedSigmaValues;
            storedSigmaValues = tempSigmaValues;
        }

        double[] tempRhoValues = rhoValues;
        rhoValues = storedRhoValues;
        storedRhoValues = tempRhoValues;

        double[] tempSigmaMatrix = sigmaMatrix;
        sigmaMatrix = storedSigmaMatrix;
        storedSigmaMatrix = tempSigmaMatrix;

        double[] tempRhoMatrix = rhoMatrix;
        rhoMatrix = storedRhoMatrix;
        storedRhoMatrix = tempRhoMatrix;

        if(inverseMatrix) {
            double[] tempInverseRhoMatrix = inverseRhoMatrix;
            inverseRhoMatrix = storedInverseRhoMatrix;
            storedInverseRhoMatrix = tempInverseRhoMatrix;
        }

        if(coEstimate) {
            double[] tempCovariance = covarianceValues;
            covarianceValues = storedCovariance;
            storedCovariance = tempCovariance;
        }
    }

    /*
     * These methods can be static
     */
    private void populateSigmaMatrixCo(double variance, double[] covariance, double[] sigma, double [] transSigma, int nTraits, double[] traitRateMatrix){
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(sigma, i, i, variance, nTraits);
            for (int j = i + 1; j < nTraits; j++) {
                MatrixUtilsContra.setMatrixEntry(sigma,i, j, covariance[k], nTraits);
                k++;
            }
        }
        MatrixUtilsContra.matrixTranspose(sigma, nTraits, transSigma);

        MatrixUtilsContra.matrixMultiply(sigma, transSigma, nTraits, nTraits, traitRateMatrix);
    }

    private void populateSigmaMatrixCo(double[] variance, double[] covariance, double[] sigma, double [] transSigma, int nTraits, double[] traitRateMatrix){
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(sigma, i, i, variance[i], nTraits);
            for (int j = i + 1; j < nTraits; j++) {
                MatrixUtilsContra.setMatrixEntry(sigma,i, j, covariance[k], nTraits);
                k++;
            }
        }
        MatrixUtilsContra.matrixTranspose(sigma, nTraits, transSigma);

        MatrixUtilsContra.matrixMultiply(sigma, transSigma, nTraits, nTraits, traitRateMatrix);
    }

    private void populateSigmaMatrixDirectlyCo(double sigmasq, double[] rho, int nTraits, double[] traitRateMatrix) {
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(traitRateMatrix, i, i, sigmasq, nTraits);
            for (int j = i + 1; j < nTraits; j++) {
                double cov = rho[k];
                MatrixUtilsContra.setMatrixEntry(traitRateMatrix,i, j, cov, nTraits);
                MatrixUtilsContra.setMatrixEntry(traitRateMatrix,j, i, cov, nTraits);
                k++;
            }
        }
    }

    private void populateSigmaMatrixDirectlyCo(double[] sigmasq, double[] rho, int nTraits, double[] traitRateMatrix) {
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(traitRateMatrix, i, i, sigmasq[i], nTraits);
            for (int j = i + 1; j < nTraits; j++) {
                double cov = rho[k];
                MatrixUtilsContra.setMatrixEntry(traitRateMatrix,i, j, cov, nTraits);
                MatrixUtilsContra.setMatrixEntry(traitRateMatrix,j, i, cov, nTraits);
                k++;
            }
        }
    }

}
