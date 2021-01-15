package contraband.math;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import org.apache.commons.math3.util.FastMath;

public class LiabilityNodeMath extends NodeMath{
    final public Input<RealParameter> inverseMatrixInput = new Input<>("inverseMatrix", "Inverse of trait correlation matrix when using Gibbs sampling.");

    private double[] inverseMatrix;
    private double[] inverseTraitRateMatrix;
    private double[] traitRateMatrix;

    // for doing LUDecompostion
    private double[] lu;
    private int[] pivot;
    private boolean[] evenSingular;
    private double [] identityMatrix;

    private int nTraits;
    private int matDim;

    private double [] storedInvTraitRateMatrix;
    private double [] storedTraitRateMatrix;
    private double[] storedInverseMatrix;


    @Override
    public void initAndValidate() {
        super.initAndValidate();
        nTraits = getNTraits();
        matDim = nTraits * nTraits;

        inverseMatrix = new double[matDim];
        inverseTraitRateMatrix = new double[matDim];
        traitRateMatrix = new double[matDim];

        storedInvTraitRateMatrix = new double[matDim];
        storedTraitRateMatrix = new double[matDim];
        storedInverseMatrix = new double[matDim];

        lu = new double [matDim];
        pivot = new int[nTraits];
        evenSingular = new boolean[2];
        identityMatrix = new double[matDim];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, nTraits);
        }

        if(inverseMatrixInput.get().getDimension() != matDim) {
            inverseMatrixInput.get().setDimension(matDim);
            // initiate an identity matrix
            for (int i = 0; i < matDim; i++){
                if(i % (nTraits + 1) == 0) {
                    inverseMatrixInput.get().setValue(i, 1.0);
                } else {
                    inverseMatrixInput.get().setValue(i, 0.0);
                }
            }
        } else {
            // get from input
            inverseMatrix = inverseMatrixInput.get().getDoubleValues();
        }

    }

    @Override
    public void populateTraitRateMatrix() {
        // since we sampling the inverse matrix
        // we only need to get the inverse matrix
        inverseMatrix = inverseMatrixInput.get().getDoubleValues();

        // then populate the inverse matrix of trait rate matrix
        double sigmaSq = getSigmaValue();
        MatrixUtilsContra.vectorMapMultiply(inverseMatrix, 1/sigmaSq, inverseTraitRateMatrix);
        setTraitRateMatrixInverse(inverseTraitRateMatrix);

        LUDecompositionForArray.ArrayLUDecomposition(inverseTraitRateMatrix, lu, pivot, evenSingular, nTraits);
        // to get the original trait rate matrix because it will be used when calculating r parameter
       LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evenSingular[1], nTraits, traitRateMatrix);
       setTraitRateMatrix(traitRateMatrix);
    }

    @Override
    public void performMatrixOperations() {
        // perform LUDecomposition and get the determinant of the inverse matrix
        operateOnInvTraitRateMatrix();
        double detInv = getTraitRateMatrixInverseDeterminant();

        // the determinant of original trait rate matrix is the reciprocal of the determinant of the inverse matrix
        double logDet = FastMath.log(1/detInv);

        setTraitRateMatrixDeterminant(logDet);
    }

    /*
     * To check if some parameters are dirty
     * for store and restor
     */
    @Override
    public boolean updateParameters() {
        boolean updateInverseMatrix = false;
        boolean others = super.updateParameters();
        if(inverseMatrixInput.isDirty()){
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
        System.arraycopy(inverseTraitRateMatrix, 0, storedInvTraitRateMatrix, 0, matDim);
        System.arraycopy(inverseMatrix, 0, storedInverseMatrix, 0, matDim);
    }

    @Override
    public void restore() {
        super.restore();

        double[] tempTraitRateMatrix = traitRateMatrix;
        traitRateMatrix = storedTraitRateMatrix;
        storedTraitRateMatrix = tempTraitRateMatrix;

        double[] tempInvTraitRateMatrix = inverseTraitRateMatrix;
        inverseTraitRateMatrix = storedInvTraitRateMatrix;
        storedInvTraitRateMatrix = tempInvTraitRateMatrix;

        double[] tempInverseMatrix = inverseMatrix;
        inverseMatrix = storedInverseMatrix;
        storedInverseMatrix = tempInverseMatrix;
    }
}
