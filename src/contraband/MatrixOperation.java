package contraband;

import org.apache.commons.math3.linear.RealMatrix;

public class MatrixOperation {

    /*
     * This method do the following calculation:
     * resRM <- inRM + rmToAdd
     * Note that: resRM is initialized and will be rewritten in the method.
     */
    public static RealMatrix matrixAdd(RealMatrix inRM, RealMatrix rmToAdd, RealMatrix resRM) {
        for (int i = 0; i < inRM.getRowDimension(); i ++) {
            for (int j = 0; j < inRM.getColumnDimension(); j ++) {
                resRM.setEntry(i, j, inRM.getEntry(i, j) + rmToAdd.getEntry(i, j));
            }
        }

        return resRM;
    }

    /*
     * aMat <- aMat - bMat
     */
    public static void matrixSubtract(RealMatrix aMat, RealMatrix bMat) {
        for (int i = 0; i < aMat.getRowDimension(); i ++) {
            for (int j = 0; j < aMat.getColumnDimension(); j ++) {
                aMat.setEntry(i, j, aMat.getEntry(i, j) - bMat.getEntry(i, j));
            }
        }
    }

    /*
     * aMat <- scalar * aMat
     */
    public static void matrixScalarMultiply(RealMatrix aMat, double scalar) {
        for (int i = 0; i < aMat.getRowDimension(); i ++) {
            for (int j = 0; j < aMat.getColumnDimension(); j ++) {
                aMat.setEntry(i, j, aMat.getEntry(i, j) * scalar);
            }
        }
    }

    /*
     * aMat <- aMat + scalar * bMat
     */
    public static void matrixScalarAdd(RealMatrix aMat, RealMatrix bMat, double scalar) {
        for (int i = 0; i < aMat.getRowDimension(); i ++) {
            for (int j = 0; j < aMat.getColumnDimension(); j ++) {
                aMat.setEntry(i, j, aMat.getEntry(i, j) + scalar * bMat.getEntry(i, j));
            }
        }
    }


    /*
     * aMat <- Mat.transpose
     */
    public static void matrixTranspose(RealMatrix aMat) {
        for (int i = 0; i < aMat.getRowDimension(); i ++) {
            for (int j = 0; j < aMat.getColumnDimension(); j ++) {
                aMat.setEntry(i, j, aMat.getEntry(j, i));
            }
        }
    }


}
