package contraband.math;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import org.apache.commons.math3.util.FastMath;

public class LiabilityNodeMath extends NodeMath{
    final public Input<RealParameter> inverseMatrixInput = new Input<>("inverseMatrix", "Inverse of trait correlation matrix when using Gibbs sampling.");

    private double[] inverseMatrix;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        int nTraits = getNTraits();
        int matDim = nTraits * nTraits;

        inverseMatrix = new double[matDim];
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
        setTraitRateMatrixInverse(inverseMatrix);
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


}
