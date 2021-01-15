package contraband.math;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import org.apache.commons.math3.util.FastMath;

public class LiabilityNodeMath extends NodeMath{
    final public Input<RealParameter> inverseRhoInput = new Input<>("inverseMatrix", "Inverse of trait correlation matrix when using Gibbs sampling.");

    private double[] inverseRho;
    private double[] inverseTraitRateMatrix;
    private double[] traitRateMatrix;
    private double[] rhoMatrix;

    // for doing LUDecompostion
    private double[] lu;
    private int[] pivot;
    private boolean[] evenSingular;
    private double [] identityMatrix;

    private int nTraits;
    private int matDim;
    private double detInvRhoMatrix;
    private double detRhoMatrix;

    private double [] storedInvTraitRateMatrix;
    private double [] storedRhoMatrix;
    private double[] storedInverseMatrix;
    private double [] storedTraitRateMatrix;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        nTraits = getNTraits();
        matDim = nTraits * nTraits;

        inverseRho = new double[matDim];
        inverseTraitRateMatrix = new double[matDim];
        traitRateMatrix = new double[matDim];
        rhoMatrix = new double[matDim];

        storedInvTraitRateMatrix = new double[matDim];
        storedRhoMatrix = new double[matDim];
        storedInverseMatrix = new double[matDim];
        storedTraitRateMatrix = new double[matDim];
        lu = new double [matDim];
        pivot = new int[nTraits];
        evenSingular = new boolean[2];
        identityMatrix = new double[matDim];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, nTraits);
        }

        if(inverseRhoInput.get().getDimension() != matDim) {
            inverseRhoInput.get().setDimension(matDim);
            // initiate an identity matrix
            for (int i = 0; i < matDim; i++){
                if(i % (nTraits + 1) == 0) {
                    inverseRhoInput.get().setValue(i, 1.0);
                } else {
                    inverseRhoInput.get().setValue(i, 0.0);
                }
            }
        } else {
            // get from input
            inverseRho = inverseRhoInput.get().getDoubleValues();
        }

    }

    @Override
    public void populateTraitRateMatrix() {
        // since we sampling the inverse matrix
        // we only need to get the inverse matrix
        inverseRho = inverseRhoInput.get().getDoubleValues();

        boolean singularMatrix = false;
        // to get the original rho matrix because it will be used when calculating r parameter
        LUDecompositionForArray.ArrayLUDecomposition(inverseRho, lu, pivot, evenSingular, nTraits);
        detInvRhoMatrix = LUDecompositionForArray.getDeterminant(lu, nTraits, evenSingular);
         try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evenSingular[1], nTraits, rhoMatrix);
        } catch (RuntimeException e) {
            singularMatrix = true;
        }
        setSingularMatrix(singularMatrix);

        if(isOneRateOnly()) {
            double sigmaSq = getSigmaValue();
            MatrixUtilsContra.vectorMapMultiply(rhoMatrix, sigmaSq, traitRateMatrix);
        } else {
            double[] sigmaSqs = getSigmaValues();
            for(int i = 0; i < nTraits; i++){
                for(int j = 0; j < nTraits; j++){
                   traitRateMatrix[i * nTraits + j] = FastMath.sqrt(sigmaSqs[i]) * FastMath.sqrt(sigmaSqs[j]) * rhoMatrix[i * nTraits + j];
                }
            }
        }
        setTraitRateMatrix(traitRateMatrix);

    }

    @Override
    public void performMatrixOperations() {
        if(isOneRateOnly()){
            // populate the inverse matrix of trait rate matrix
            double sigmaSq = getSigmaValue();

            MatrixUtilsContra.vectorMapMultiply(inverseRho, 1/ sigmaSq, inverseTraitRateMatrix);
            setInvTraitRateMatrix(inverseTraitRateMatrix);

            // the determinant of original trait rate matrix is the reciprocal of the determinant of the inverse matrix
            detRhoMatrix = (1.0 / sigmaSq) * detInvRhoMatrix;

            double detInvTraitRateMat = detRhoMatrix * FastMath.pow(1/sigmaSq, nTraits);
            setInvTraitRateMatrixDeterminant(detInvTraitRateMat);

            double logDetTraitRateMat = FastMath.log(detRhoMatrix) + nTraits * FastMath.log(sigmaSq);
            setTraitRateMatrixDeterminant(logDetTraitRateMat);
        } else {
            // for multiple trait evolutionary rates
            // we will need to LUDecomposition and get the determinant of the matrix
            super.performMatrixOperations();
        }
    }

    /*
     * To check if some parameters are dirty
     * for store and restor
     */
    @Override
    public boolean updateParameters() {
        boolean updateInverseMatrix = false;
        boolean others = super.updateParameters();
        if(inverseRhoInput.isDirty()){
            updateInverseMatrix = true;
        }
        return others || updateInverseMatrix;
    }

    /*
     * Since we sample the inverse matrix, we don't need to deal with Rho matrix any more.
     */
    @Override
    protected void initiateRhoInput(){ }


    @Override
    public void store() {
        super.store();
        System.arraycopy(traitRateMatrix, 0, storedTraitRateMatrix, 0, matDim);
        System.arraycopy(rhoMatrix, 0, storedRhoMatrix, 0, matDim);
        System.arraycopy(inverseTraitRateMatrix, 0, storedInvTraitRateMatrix, 0, matDim);
        System.arraycopy(inverseRho, 0, storedInverseMatrix, 0, matDim);
    }

    @Override
    public void restore() {
        super.restore();

        double[] tempRhoMatrix = rhoMatrix;
        rhoMatrix = storedRhoMatrix;
        storedRhoMatrix = tempRhoMatrix;

        double[] tempInvTraitRateMatrix = inverseTraitRateMatrix;
        inverseTraitRateMatrix = storedInvTraitRateMatrix;
        storedInvTraitRateMatrix = tempInvTraitRateMatrix;

        double[] tempInverseMatrix = inverseRho;
        inverseRho = storedInverseMatrix;
        storedInverseMatrix = tempInverseMatrix;

        double[] tempTraitRateMatrix = traitRateMatrix;
        traitRateMatrix = storedTraitRateMatrix;
        storedTraitRateMatrix = tempTraitRateMatrix;
    }
}
