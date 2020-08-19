package contraband.math;

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
    
    public static void elementWiseProduct(double[][] mat1, double[][] mat2, double[][] res) {
        int rows = mat1.length;
        int cols = mat1[0].length;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                res[i][j] = mat1[i][j] * mat2[i][j];
            }
        }
    }

    // For computing the EB covariance matrix and OU multivariate case
    public static void kronecker(double[][] mat1, double[][] mat2, double[][] res) {
        int r1 = mat1.length;
        int c1 = mat1[0].length;
        int r2 = mat2.length;
        int c2 = mat2[0].length;

        for (int i1 = 0; i1 < r1; i1++) {
            for (int j1 = 0; j1 < c1; j1++) {
                for (int i2 = 0; i2 < r2; i2++) {
                    for (int j2 = 0; j2 < c2; j2++) {
                        res[i1 * r2 + i2][j1 * c2 + j2] = mat1[i1][j1] * mat2[i2][j2];
                    }
                }
            }
        }
    }

    public static void scalarByRealMatrix(RealMatrix aRealMat, double scalar) {
        int rows = aRealMat.getRowDimension();
        int cols = aRealMat.getColumnDimension();

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                aRealMat.multiplyEntry(i, j, scalar);
            }
        }
    }

    public static int[] forLoopRange(int from, int limit) {
        int[] numbers = new int[limit] ;
        for (int i = 0; i < limit; i++) {
            numbers[i] = from;
            from++;
        }

        return numbers;
    }

    public static void setBlocksInMatrix(double[][] blockMat, int m, int rowIndex, int colIndex, double[][] Result) {

        int i = 0, j = 0;
        int[] rowElements = forLoopRange(rowIndex * m, m); System.out.println(" Numero de filas filas " + rowElements.length);
        int[] colElements = forLoopRange(colIndex * m, m); System.out.println(" Numero de columnas " + colElements.length);

        for(int row:rowElements) {

            for(int col:colElements) {

                System.out.println(col);
                Result[row][col] = 	blockMat[i][j]; System.out.println("row = " + row + "," + "col = " + col);

                j++;
            }

            j = 0;
            i++;
        }
    }
}
