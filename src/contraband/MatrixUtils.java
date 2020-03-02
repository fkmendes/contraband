package contraband;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class MatrixUtils {

    /*
     * This method adds matrix inRM to rmToAdd.
     *
     * Note that this method returns a RealMatrix
     * resRM, but this matrix was declared and initialized
     * elsewhere -- meaning this method has a side-effect.
     */
    public static RealMatrix matrixAdd(RealMatrix inRM, RealMatrix rmToAdd, RealMatrix resRM) {
        if (inRM.getColumnDimension() == rmToAdd.getColumnDimension()
            && inRM.getRowDimension() == rmToAdd.getRowDimension()
            && inRM.getColumnDimension() == resRM.getColumnDimension()
            && inRM.getRowDimension() == resRM.getRowDimension()) {
            for (int i = 0; i < inRM.getRowDimension(); i++) {
                for (int j = 0; j < inRM.getColumnDimension(); j++) {
                    resRM.setEntry(i, j, inRM.getEntry(i, j) + rmToAdd.getEntry(i, j));
                }
            }
        } else {
            throw new IllegalArgumentException("Dimension does not match in matrix add!");
        }
        return resRM;
    }

    /*
     * This method subtracts matrix rmToSubtract from inRM.
     *
     * Note that this method returns a RealMatrix
     * resRM, but this matrix was declared and initialized
     * elsewhere -- meaning this method has a side-effect.
     */
    public static RealMatrix matrixSubtract(RealMatrix inRM, RealMatrix rmToSubtract, RealMatrix resRM) {
        if (inRM.getColumnDimension() == rmToSubtract.getColumnDimension()
                && inRM.getRowDimension() == rmToSubtract.getRowDimension()
                && inRM.getColumnDimension() == resRM.getColumnDimension()
                && inRM.getRowDimension() == resRM.getRowDimension()) {
            for (int i = 0; i < inRM.getRowDimension(); i++) {
                for (int j = 0; j < inRM.getColumnDimension(); j++) {
                    resRM.setEntry(i, j, inRM.getEntry(i, j) - rmToSubtract.getEntry(i, j));
                }
            }
        } else {
            throw new IllegalArgumentException("Dimension does not match in matrix subtract!");
        }
        return resRM;
    }

    /*
     *
     */
    public static RealMatrix matrixMultiply(RealMatrix inRM, RealMatrix rmToMultiply, RealMatrix resRM) {
        if (inRM.getColumnDimension() == rmToMultiply.getRowDimension()
            && inRM.getRowDimension() == resRM.getRowDimension()
            && rmToMultiply.getColumnDimension() == resRM.getColumnDimension()) {
            for (int i = 0; i < resRM.getRowDimension(); i ++) {

                for (int j = 0; j < resRM.getColumnDimension(); j ++) {
                    double sumTemp = 0.0;

                    for (int k = 0; k < rmToMultiply.getRowDimension(); k++) {
                        sumTemp += inRM.getEntry(i, k) * rmToMultiply.getEntry(k, j);
                    }

                    resRM.setEntry(i, j, sumTemp);
                }
            }
        } else {
            throw new IllegalArgumentException("Dimension does not match in matrix multiply!");
        }
        return resRM;
    }

    /*
     * This method multiplies matrix inRM by a double scalar.
     *
     * Note that this method returns a RealMatrix
     * resRM, but this matrix was declared and initialized
     * elsewhere -- meaning this method has a side-effect.
     */
    public static RealMatrix matrixScalarMultiply(RealMatrix inRM, double scalar, RealMatrix resRM) {
        if (inRM.getColumnDimension() == resRM.getColumnDimension()
            && inRM.getRowDimension() == resRM.getRowDimension()) {
            for (int i = 0; i < inRM.getRowDimension(); i++) {
                for (int j = 0; j < inRM.getColumnDimension(); j++) {
                    resRM.setEntry(i, j, inRM.getEntry(i, j) * scalar);
                }
            }
        } else {
            throw new IllegalArgumentException("Dimension does not match in matrix scalar multiply!");
        }
        return resRM;
    }

    /*
     * This method adds a double scalar to matrix inRM.
     *
     * Note that this method returns a RealMatrix
     * resRM, but this matrix was declared and initialized
     * elsewhere -- meaning this method has a side-effect.
     */
    public static RealMatrix matrixScalarAdd(RealMatrix inRM, double scalar, RealMatrix resRM) {
        if (inRM.getColumnDimension() == resRM.getColumnDimension()
                && inRM.getRowDimension() == resRM.getRowDimension()) {
            for (int i = 0; i < inRM.getRowDimension(); i++) {
                for (int j = 0; j < inRM.getColumnDimension(); j++) {
                    resRM.setEntry(i, j, inRM.getEntry(i, j) + scalar);
                }
            }
        } else {
            throw new IllegalArgumentException("Dimension does not match in matrix scalar add!");
        }
        return resRM;
    }


    /*
     * This method transposes matrix inRM.
     *
     * Note that this method returns a RealMatrix
     * resRM, but this matrix was declared and initialized
     * elsewhere -- meaning this method has a side-effect.
     */
    public static RealMatrix matrixTranspose(RealMatrix inRM, RealMatrix resRM) {
        if (inRM.getColumnDimension() == resRM.getRowDimension()
            && inRM.getRowDimension() == resRM.getColumnDimension()) {
            for (int i = 0; i < inRM.getRowDimension(); i++) {
                for (int j = 0; j < inRM.getColumnDimension(); j++) {
                    resRM.setEntry(j, i, inRM.getEntry(i, j));
                }
            }
        } else {
            throw new IllegalArgumentException("Dimension does not match in matrix transpose!");
        }
        return resRM;
    }

    /*
     * This method calculates vector dotProduct and matrix preMultiply, i.e.
     * inVec.transpose * inRm * inVec,
     * and returns a double value.
     */
    public static double vecDotMatrixPreMultiply(RealVector inVec, RealMatrix inRM) {
        double res = 0.0;
        for (int i = 0; i < inRM.getColumnDimension(); i++) {
            double sum = 0.0;
            for (int j = 0; j < inRM.getRowDimension(); j++) {
                sum = sum + inVec.getEntry(j) * inRM.getEntry(i, j);
            }
            res = res + sum * inVec.getEntry(i);
        }
        return res;
    }

    /*
     * This method calculates matrix preMultiply, i.e.
     * resVec <- inVec * inRM
     *
     * Note that this method not only populates resVec,
     * but also returns resVec.
     */
    public static RealVector matrixPremultiply(RealVector inVec, RealMatrix inRM, RealVector resVec) {
        for (int i = 0; i < inRM.getColumnDimension(); i++) {
            double sum = 0.0;
            for (int j = 0; j < inRM.getRowDimension(); j++) {
                sum = sum + inVec.getEntry(j) * inRM.getEntry(i, j);
            }
            resVec.setEntry(i, sum);
        }
        return resVec;
    }

    /*
     * This method calculates vector map multiply, i.e.
     * resVec <- inVec * scalar
     *
     * Note that this not only populates resVec,
     * but also returns resVec.
     */
    public static RealVector vectorMapmultiply(RealVector inVec, double scalar, RealVector resVec) {
        for (int i = 0; i < inVec.getDimension(); i++) {
            resVec.setEntry(i, inVec.getEntry(i) * scalar);
        }
        return resVec;
    }

    /*
     * This method adds inVec to vecToAdd, i.e.
     * resVec <- inVec + vecToAdd
     */
    public static RealVector vectorAdd (RealVector inVec, RealVector vecToAdd, RealVector resVec) {
        for (int i = 0; i < inVec.getDimension(); i++) {
            resVec.setEntry(i, inVec.getEntry(i) + vecToAdd.getEntry(i));
        }
        return resVec;
    }

}
