package contraband.math;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import contraband.utils.NodeMathUtils;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.util.FastMath;
import outercore.parameter.KeyRealParameter;

public class NodeMath extends CalculationNode {
    final public Input<KeyRealParameter> traitsValuesInput = new Input<>("traits","Trait values at tips.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Evolutionary rates of traits. Diagonal elements in rate matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> rhoInput = new Input<>("correlation", "Correlations of traits. Off-diagonal elements in rate matrix.");
    final public Input<Boolean> oneRateOnlyInput = new Input<>("oneRateOnly", "TRUE, if all traits share one evolutionary rate.", false);
    final public Input<Boolean> useUpperMatrixInput = new Input<>("upperMatrix", "TRUE, if sigmasq and correlations are from upper matrix", false);
    final public Input<RealParameter> rootValuesInput = new Input<>("rootValues", "Trait values at the root.");
    final public Input<Boolean> useShrinkageInput = new Input<>("shrinkage", "TRUE, if shrinkage method is used to estimate the trait correlations.", false);
    final public Input<RealParameter> covarianceInput = new Input<>("covariance", "cov_ij = sigma_i * sigma_j * rho_ij.");

    private Integer nTraits;
    private Integer nSpecies;
    private Boolean oneRateOnly;
    private Boolean useUpperMatrix;
    private Boolean sampleRoot;
    private Boolean coEstimate;

    private int nodeNr;
    private int vecArrayDim;
    private int matDim;

    // BM parameters
    private double sigmaValue;
    private double[] sigmaValues;
    private double[] rhoValues;


    // for doing LUDecompostion
    private double[] lu;
    private int[] pivot;
    private boolean[] evensingular;
    private double [] identityMatrix;
    private boolean singularMatrix;

    // PCM model parameters
    private double [] aArray;
    private double [] cArray;
    private double [] eArray;
    private double [] fArray;
    private double [] lArray;
    private double [] mVecArray;
    private double [] rArray;
    private double [] traitRateMatrix;
    private double detTraitRateMat;
    private double [] invTraitRateMatrix;
    private double detInvTraitRateMat;

    // temporary parameters
    private double [] mVecInit;
    private double[] traitsVec;
    private double [] rateMatRow;
    private double [] upperMatrix;
    private double [] transUpperMatrix;
    private double [] mVec;

    // parameters at the root
    private double [] rootValuesVec;

    // this is for trees with sampled ancestors
    private double likForSA;
    private double[] traitVecForSA;

    // variance and expectation for BM
    private double [] nodeVariance;
    private double [] nodeExpectation;
    private double [] expect1;
    private double [] expect2;
    private double [] expectp;

    private double storedDetTraitRateMat;
    private double [] storedInvTraitRateMatrix;
    private double storedDetInvTraitRateMat;
    private double [] storedRootValuesVec;
    private double storedSigmaValue;
    private double[] storedSigmaValues;
    private double[] storedRhoValues;
    private double [] storedTraitRateMatrix;

    // shrinkage
    private Boolean useShrinkage;
    private RealMatrix unbiasedRhoRM;
    private RealMatrix identityRM;
    private RealMatrix traitRateRM;
    private RealMatrix shrinkageRhoRM;
    private double detRhoMatrix;
    private double detInvRhoMatrix;
    private double[] transformedTraitValues;
    private double[] storedTransformedTraitValues;


    @Override
    public void initAndValidate() {
        // collect trait information
        KeyRealParameter traitsValues = traitsValuesInput.get();
        nTraits = traitsValues.getMinorDimension1();
        nSpecies = traitsValues.getMinorDimension2();

        // TRUE, if sigmasq and rho are in variance-covariance parameterization.
        if (covarianceInput.get() == null) {
            coEstimate = false;
        } else {
            coEstimate = true;
            Log.info("NodeMath::Variance-covariance parameterization is used.");

            // check dimension
            if(covarianceInput.get().getDimension() != (nTraits * (nTraits - 1) / 2)) {
                Log.warning.println("NodeMath::Setting dimension of covariance to " + (nTraits * (nTraits - 1) / 2));
                covarianceInput.get().setDimension((nTraits * (nTraits - 1) / 2));
            }
            rhoValues = new double[(nTraits * (nTraits - 1) / 2)];
            storedRhoValues = new double[(nTraits * (nTraits - 1) / 2)];
            rhoValues = covarianceInput.get().getDoubleValues();
        }


        // initialize real matrix
        identityRM = MatrixUtils.createRealIdentityMatrix(nTraits);
        traitRateRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        shrinkageRhoRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        transformedTraitValues = new double[nSpecies * nTraits];
        storedTransformedTraitValues = new double[nSpecies * nTraits];

        // default is false, i.e. each trait has its own rate
        oneRateOnly = oneRateOnlyInput.get();

        // default is false, i.e. not doing t(Sigma) * Sigma
        useUpperMatrix = useUpperMatrixInput.get();

        // dimensions
        nodeNr = 2 * nSpecies - 1;
        matDim = nTraits * nTraits;
        vecArrayDim = nTraits * nodeNr;

        // initialize matrix and vectors
        initializeAbCdEfArray();
        initializeLmrArray();
        initializeTraitValueVec();

        traitRateMatrix = new double[matDim];
        invTraitRateMatrix = new double[matDim];
        storedTraitRateMatrix = new double[matDim];
        storedInvTraitRateMatrix = new double[matDim];

        rateMatRow = new double [nTraits];
        if(useUpperMatrix) {
            upperMatrix = new double[matDim];
            transUpperMatrix = new double[matDim];
        }
        // each node has its variance = branch rate * branch length
        nodeVariance = new double [nodeNr];
        // each node has its expectation = weight average of trait values from descendants
        nodeExpectation = new double [nTraits * nodeNr];
        // initialize expectation
        expect1 = new double [nTraits];
        expect2 = new double [nTraits];
        expectp = new double [nTraits];

        initializeLUDecomposition();

        likForSA = 0.0;
        traitVecForSA = new double [nTraits];

        useShrinkage = useShrinkageInput.get();
        if (!useShrinkage && !coEstimate) {
            if (rhoInput.get() == null) {
                throw new RuntimeException("NodeMath::If shrinkage method is not used, either correlation or covariance is required.");
            } else {
                Log.warning.println("NodeMath::Setting dimension of rho values to " + (nTraits * (nTraits - 1) / 2));
                rhoInput.get().setDimension((nTraits * (nTraits - 1) / 2));
                rhoValues = new double[(nTraits * (nTraits - 1) / 2)];
                storedRhoValues = new double[(nTraits * (nTraits - 1) / 2)];
                rhoValues = rhoInput.get().getDoubleValues();
            }
        }

        // check rate input
        if (oneRateOnly) {
            Log.info.println("NodeMath::1 rate is assigned to " + nTraits + " traits.");
            sigmasqInput.get().setDimension(1);
            sigmaValue = sigmasqInput.get().getValue();
            populateTraitRateMatrixForOneRate();
        } else {
            Log.warning.println("NodeMath::Setting dimension of trait rate to " + nTraits + ".");
            sigmasqInput.get().setDimension(nTraits);
            sigmaValues = new double[nTraits];
            storedSigmaValues = new double[nTraits];
            sigmaValues = sigmasqInput.get().getDoubleValues();
            populateTraitRateMatrixForMultipleRates();
        }

        // check root values dimension
        if(rootValuesInput.get() != null) {
            if (rootValuesInput.get().getDimension() != nTraits) {
                Log.warning.println("NodeMath::Setting dimension of root values to " + nTraits);
                rootValuesInput.get().setDimension(nTraits);
            }
            rootValuesVec = rootValuesInput.get().getDoubleValues();
            sampleRoot = true;
        } else {
            Log.warning.println("NodeMath::Estimating root state by MLE.");
            sampleRoot = false;
        }
    }

    //getters
    public boolean isSingularMatrix () { return singularMatrix; }

    public double getAForNode (int nodeIdx) { return aArray[nodeIdx]; }

    public double getCForNode (int nodeIdx) { return cArray[nodeIdx]; }

    public double getEForNode (int nodeIdx) { return eArray[nodeIdx]; }

    public double getfForNode (int nodeIdx) { return fArray[nodeIdx]; }

    public double getLForNode (int nodeIdx) { return lArray[nodeIdx]; }

    public double getRForNode (int nodeIdx) { return rArray[nodeIdx]; }

    public double[] getMVecForNode (int nodeIdx) {
        MatrixUtilsContra.getRowArray(mVecArray, nodeIdx, nTraits, mVec);
        return mVec;
    }

    public double[] getTraitRateMatrix () { return traitRateMatrix; }

    public double getTraitRateMatrixDeterminant () { return detTraitRateMat; }

    public double[] getTraitRateMatrixInverse () { return invTraitRateMatrix; }

    public double getTraitRateMatrixInverseDeterminant () { return detInvTraitRateMat; }

    public double[] getRootValuesArr () { return rootValuesVec; }

    public double getLikelihoodForSampledAncestors () { return likForSA; }

    public double[] getTempVec () { return mVec; }

    public double[] getInitMVec () { return mVecInit; }

    public double[] getTraitsVec () {return traitsVec; }

    public double[] getSampledAncestorTraitsVec () { return traitVecForSA; }

    public double[] getRateMatrixRow () { return rateMatRow; }

    public double getVarianceForNode (int nodeIdx) { return nodeVariance[nodeIdx]; }

    public double [] getTransformedTraitValues () { return transformedTraitValues; }

    public boolean useShrinkage () { return useShrinkage; }

    //setters
    public void setAForNode (int nodeIdx, double value) { aArray[nodeIdx] = value; }

    public void setCForNode (int nodeIdx, double value) { cArray[nodeIdx] = value; }

    public void setEForNode (int nodeIdx, double value) { eArray[nodeIdx] = value; }

    public void setfForNode (int nodeIdx, double value) { fArray[nodeIdx] = value; }

    public void setLForNode (int nodeIdx, double value) { lArray[nodeIdx] = value; }

    public void setRForNode (int nodeIdx, double value) { rArray[nodeIdx] = value; }

    public void setMVecForNode (int nodeIdx, double[] value) {
        MatrixUtilsContra.setArray(mVecArray, value, nodeIdx, nTraits);
    }

    public void setLikelihoodForSampledAncestors (double value) { likForSA = value; }

    public void setTraitsVecForTip (double[] traitValues, int tipIdx) {
        MatrixUtilsContra.getRowArray(traitValues, tipIdx, nTraits, traitsVec);
    }

    public void setTraitsVecForSampledAncestor (double[] traitValues, int nodeIdx) {
        MatrixUtilsContra.getRowArray(traitValues, nodeIdx, nTraits, traitVecForSA);
    }

    public void setExpectationForTip (int nodeIdx) {
        MatrixUtilsContra.setArray(nodeExpectation, traitsVec, nodeIdx, nTraits);
    }

    public void setVarianceForTip (int nodeIdx, double value) { nodeVariance[nodeIdx] = value; }

    public void setVarianceForParent (int parentIdx, double branchLength, int child1Idx, int child2Idx) {
        double vc1 = nodeVariance[child1Idx];
        double vc2 = nodeVariance[child2Idx];
        nodeVariance[parentIdx] = branchLength + ((vc1 * vc2)/(vc1 + vc2));
    }

    public void setExpectationForParent (int parentIdx, int child1Idx, int child2Idx) {
        double vc1 = nodeVariance[child1Idx];
        double vc2 = nodeVariance[child2Idx];
        MatrixUtilsContra.getRowArray(nodeExpectation, child1Idx, nTraits, expect1);
        MatrixUtilsContra.getRowArray(nodeExpectation, child2Idx, nTraits, expect2);
        MatrixUtilsContra.vectorMapMultiply(expect1, vc2/(vc1+vc2), expect1);
        MatrixUtilsContra.vectorMapMultiply(expect2, vc1/(vc1+vc2), expect2);
        MatrixUtilsContra.vectorAdd(expect1, expect2, expectp);
        MatrixUtilsContra.setArray(nodeExpectation, expectp, parentIdx, nTraits);
    }

    public void setNTraits (int value) { nTraits = value; }

    public void setNSpecies (int value) { nSpecies = value; }

    private void initializeAbCdEfArray() {
        // A C E f are single double values for each node
        aArray = new double[nodeNr];
        cArray = new double[nodeNr];
        eArray = new double[nodeNr];
        fArray = new double[nodeNr];
    }

    private void initializeLmrArray() {
        // each node has a single double value L and r
        // and a vector m
        lArray = new double[nodeNr];
        mVecArray = new double[vecArrayDim];
        rArray = new double[nodeNr];

        mVec = new double [nTraits];
        mVecInit = new double [nTraits];
    }


    private void initializeTraitValueVec() {
        rootValuesVec = new double [nTraits];
        storedRootValuesVec = new double [nTraits];
        traitsVec = new double [nTraits];
    }

    private void initializeLUDecomposition() {
        lu = new double [matDim];
        pivot = new int[nTraits];
        evensingular = new boolean[2];
        identityMatrix = new double[matDim];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, nTraits);
        }
    }


    public void populateRootValuesVec(int rootIdx) {
        if(sampleRoot && rootValuesInput.isDirty()) {
            rootValuesVec = rootValuesInput.get().getDoubleValues();
        }
        if(!sampleRoot) {
            MatrixUtilsContra.getRowArray(nodeExpectation, rootIdx, nTraits, rootValuesVec);
        }
    }

    /*
     * This method performs matrix operations to get the determinant and inverse matrix.
     * (1)use shrinkage
     * -> populates the trait rate matrix by estimated correlations
     * (2)not use shrinkage
     * -> populates the trait rate matrix by rho and sigmasq inputs
     *
     * (3) matrix operations
     * -> inverse of trait rate matrix:  invTraitRateMatrx
     * -> determinant of the trait rate matrix: detTraitRateMat
     * -> determinant of the inverse trait rate matrix: detInvTraitRateMat
     */
    public void performMatrixOperations() {
            // determinant and inverse of trait rate matrix
            operateOnTraitRateMatrix();
            // determinant of inverse trait rate matrix
            operateOnInvTraitRateMatrix();
    }

    /*
     * This method implements LUDecomposition
     * and calculates determinant and inverse of trait rate matrix.
     *
     * Note: trait rate matrix is not populated by shrinkage method.
     */
    private void operateOnTraitRateMatrix(){
        singularMatrix = false;

        // LUDecomposition
        LUDecompositionForArray.ArrayLUDecomposition(traitRateMatrix, lu, pivot, evensingular, nTraits);

        // invert the traitRateMatrix
        try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evensingular[1], nTraits, invTraitRateMatrix);
        } catch (RuntimeException e) {
            singularMatrix = true;
        }
        // get the determinant of traitRateMatrix
        double det = LUDecompositionForArray.getDeterminant(lu, nTraits, evensingular);
        if (det == 0.0) {
            singularMatrix = true;
        }
        detTraitRateMat = FastMath.log(det);
    }

    /*
     * This method implements LUDecomposition
     * and calculates determinant of inverse trait rate matrix.
     *
     * Note: trait rate matrix is not populated by shrinkage method.
     */

    private void operateOnInvTraitRateMatrix() {
        // LUDecomposition
        LUDecompositionForArray.ArrayLUDecomposition(invTraitRateMatrix, lu, pivot, evensingular, nTraits);

        // get the determinant of invTraitRateMatrix
        detInvTraitRateMat = LUDecompositionForArray.getDeterminant(lu, nTraits, evensingular);
    }

    /*
     * This method populates trait rate matrix by trait correlations and trait evolutionary rates.
     *
     * Note: this method is applied to both using shrinkage and not using shrinkage method.
     *
     */
    public void populateTraitRateMatrix() {
        if(oneRateOnly) {
            // one rate for all trait
            // trait rate matrix is a scalar product of Rho matrix and sigmasq.
            populateTraitRateMatrixForOneRate();
        } else {
            // each trait has its own trait
            // diagonal elements = sigmasq_i
            // off-diagonal elements = sigmasq_i * sigmasq_j * rho_ij
            populateTraitRateMatrixForMultipleRates();
        }

    }

    /*
     * This method populates trait rate matrix when oneRateOnly = TRUE.
     *
     * (1) use shrinkage
     * -> trait rate matrix = sigmasq * Rho matrix
     * -> inverse trait rate matrix = (1/sigmasq) * Rho matrix
     * where Rho matrix is estimated by shrinkage method
     *
     * (2) not use shrinkage: trait rate matrix is populated by input rho and sigmasq.
     *
     */
    public void populateTraitRateMatrixForOneRate() {
        if(useShrinkage){

            traitRateRM = shrinkageRhoRM.scalarMultiply(sigmaValue);

        } else {

            if (useUpperMatrix) {
                // trait rate matrix = t(upperMatrix) * upperMatrix
                // where upperMatrix is from rho and sigmasq
                NodeMathUtils.populateTraitRateMatrix(sigmaValue, rhoValues, upperMatrix, transUpperMatrix, nTraits, traitRateMatrix);
            } else {
                // trait rate matrix_ji =  trait rate matrix_ij
                NodeMathUtils.populateTraitRateMatrixDirectly(sigmaValue, rhoValues, nTraits, traitRateMatrix);
            }
        }
    }

    /*
     * This method populates trait rate matrix when oneRateOnly = FALSE.
     *
     * (1) use shrinkage
     * -> trait rate matrix is populated by Rho matrix from shrinkage method and input sigmasq values
     * -> sqrt(sigmasq_i) * sqrt(sigmasq_j) * Rho_ij
     *
     * (2) not use shrinkage
     * -> trait rate matrix is populated by input rho and sigmasq.
     * -> sqrt(sigmasq_i) * sqrt(sigmasq_j) * rho_x
     * -> where rho_x represents the correlation between trait i and trait j
     *
     */
    public void populateTraitRateMatrixForMultipleRates() {
        if(useShrinkage) {

            for (int i = 0; i < nTraits; i++) {
                traitRateRM.setEntry(i, i, sigmaValues[i]);
                for (int j = i + 1; j < nTraits; j++) {
                    double cov = FastMath.sqrt(sigmaValues[i]) * FastMath.sqrt(sigmaValues[j]) * shrinkageRhoRM.getEntry(i, j);
                    traitRateRM.setEntry(i, j, cov);
                }
            }

        } else {

            if (useUpperMatrix) {
                // trait rate matrix = t(upperMatrix) * upperMatrix
                NodeMathUtils.populateTraitRateMatrix(sigmaValues, rhoValues, upperMatrix, transUpperMatrix, nTraits, traitRateMatrix, coEstimate);
            } else {
                // trait rate matrix_ji =  trait rate matrix_ij
                NodeMathUtils.populateTraitRateMatrixDirectly(sigmaValues, rhoValues, nTraits, traitRateMatrix);
            }

        }
    }

    /*
     * This method updates the parameters in BM model.
     *
     * In mcmc chain, if the parameters is dirty,
     * the likelihood will get the latest values.
     */
    public boolean updateParameters(){
        boolean updateRho = false;
        boolean updateSigma = false;

        // update trait correlations
        if(rhoInput.isDirty()) {
            rhoValues = rhoInput.get().getDoubleValues();
            updateRho = true;
        }
        // update trait correlations
        if(covarianceInput.isDirty()) {
            rhoValues = covarianceInput.get().getDoubleValues();
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
        return updateRho || updateSigma;
    }

    // This method calculates unbiased estimation of correlation matrix
    public void estimateCorrelations(RealMatrix traitMat){
        PearsonsCorrelation correlation = new PearsonsCorrelation(traitMat);
        unbiasedRhoRM = correlation.getCorrelationMatrix();
    }

    //This method implements the linear shrinkage method
    public void populateShrinkageEstimation(double delta) {
        // R = delta * I + (1 - delta) * R^
        shrinkageRhoRM = identityRM.scalarMultiply(delta).add(unbiasedRhoRM.scalarMultiply(1 - delta));

        // get the inverse and determinant of the shrinkage correlation matrix
        LUDecomposition corrMatLUD = new LUDecomposition(shrinkageRhoRM);
        RealMatrix invShrinkageCorrRM = corrMatLUD.getSolver().getInverse();
        detRhoMatrix = corrMatLUD.getDeterminant();

        // get the determinant of the inverse shrinkage correlation matrix
        LUDecomposition invCorrMatLUD = new LUDecomposition(invShrinkageCorrRM);
        detInvRhoMatrix= invCorrMatLUD.getDeterminant();
    }

    // This method populates determinant based on Rho matrix
    public void populateInverseTraitRateMatrx() {
        double sigmasq = sigmasqInput.get().getValue();
        detTraitRateMat = Math.log(detRhoMatrix) + nTraits * FastMath.log(sigmasq);
        detInvTraitRateMat = Math.log(detInvRhoMatrix) + nTraits * FastMath.log(1 / sigmasq);
    }

    // This method transforms original data set so that traits are independent of each other
    // X.transpose * V.inverse * X ---> Z.transpose * Z
    public void populateTransformedTraitValues(RealMatrix traitMat) {
        CholeskyDecomposition corrMatChol = new CholeskyDecomposition(traitRateRM);
        RealMatrix upperMat = corrMatChol.getLT();
        LUDecomposition upperMatLUD = new LUDecomposition(upperMat);
        RealMatrix dataTransformMat = upperMatLUD.getSolver().getInverse();

        RealMatrix transformedTraitsMat = traitMat.multiply(dataTransformMat);
        for (int j = 0; j < traitMat.getRowDimension(); j++) {
            System.arraycopy(transformedTraitsMat.getRow(j), 0, transformedTraitValues, j * nTraits, nTraits);
        }
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }
    @Override
    public void store() {
        super.store();
        storedDetTraitRateMat = detTraitRateMat;
        storedDetInvTraitRateMat = detInvTraitRateMat;

        if(oneRateOnly){
            storedSigmaValue = sigmaValue;
        } else {
            System.arraycopy(sigmaValues, 0, storedSigmaValues, 0, nTraits);
        }

        System.arraycopy(rhoValues, 0, storedRhoValues, 0, (nTraits * (nTraits - 1) / 2));
        System.arraycopy(traitRateMatrix, 0, storedTraitRateMatrix, 0, matDim);
        System.arraycopy(invTraitRateMatrix, 0, storedInvTraitRateMatrix, 0, matDim);
        System.arraycopy(rootValuesVec, 0, storedRootValuesVec, 0, nTraits);

        // shrinkage only
        System.arraycopy(transformedTraitValues, 0, storedTransformedTraitValues, 0, nSpecies * nTraits);
    }

    @Override
    public void restore() {
        super.restore();

        double tempDetTraitRateMat = detTraitRateMat;
        detTraitRateMat = storedDetTraitRateMat;
        storedDetTraitRateMat = tempDetTraitRateMat;

        double tempDetInvTraitRateMat = detInvTraitRateMat;
        detInvTraitRateMat = storedDetInvTraitRateMat;
        storedDetInvTraitRateMat = tempDetInvTraitRateMat;

        double[] tempTraitRateMatrix = traitRateMatrix;
        traitRateMatrix = storedTraitRateMatrix;
        storedTraitRateMatrix = tempTraitRateMatrix;

        double[] tempInvTraitRateMatrix = invTraitRateMatrix;
        invTraitRateMatrix = storedInvTraitRateMatrix;
        storedInvTraitRateMatrix = tempInvTraitRateMatrix;

        double[] tempRootValuesVec = rootValuesVec;
        rootValuesVec= storedRootValuesVec;
        storedRootValuesVec = tempRootValuesVec;

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

        // shrinkage only
        double[] tempTraitValuesArrayList = transformedTraitValues;
        transformedTraitValues = storedTransformedTraitValues;
        storedTransformedTraitValues = tempTraitValuesArrayList;
    }
}