package contraband.math;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.Arrays;

public class MatrixUtilsContra {

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
     * inVec.transpose * inRm * inVec, and returns a double value.
     */
    public static double tVecDotMatrixDotVec(RealVector inVec, RealMatrix inRM) {
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
     * This method calculates the dot product of
     * a row vector, a matrix, and the same vector as column vector,
     *
     * i.e., if input vector is column vector, this method
     * returns
     * t(vec) *dot* mat *dot* vec
     */
    public static double tVecDotMatrixDotVec(double[] inArr, double[] inMatArr, int nCol) {
        double res = 0.0;

        for (int i = 0; i < nCol; i++) {
            double sum = 0.0;

            for (int j = 0; j < nCol; j++) {
                sum = sum + inArr[j] * inMatArr[i * nCol + j];
            }

            res = res + sum * inArr[i];
        }

        return res;
    }

    /*
     * This method calculates matrix preMultiply, i.e.,
     * resVec <- inVec * inRM
     *
     * Populates and returns resVec
     */
    public static RealVector matrixPreMultiply(RealVector inVec, RealMatrix inRM, RealVector resVec) {
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
     * This method calculates a row vector-matrix product, i.e.,
     * resVec <- inVec * inRM
     *
     * Populates and returns resArr.
     */
    public static double[] matrixPreMultiply(double[] inArr, double[] inMatArr, int nCol, int nRow, double[] resArr) {
        for (int i = 0; i < nCol; i++) {
            double sum = 0.0;

            for (int j = 0; j < nRow; j++) {
                sum = sum + inArr[j] * inMatArr[j * nCol + i];
            }

            resArr[i] = sum;
        }

        return resArr;
    }

    /*
     * This method calculates vector map multiply, i.e.
     * resVec <- inVec * scalar
     *
     * Populates and returns resVec.
     */
    public static RealVector vectorMapMultiply(RealVector inVec, double scalar, RealVector resVec) {
        for (int i = 0; i < inVec.getDimension(); i++) {
            resVec.setEntry(i, inVec.getEntry(i) * scalar);
        }

        return resVec;
    }

    /*
     * This method calculates vector map multiply, i.e.
     * resArr <- inArr * scalar
     *
     * Populates and returns resArr.
     */
    public static double[] vectorMapMultiply(double[] inArr, double scalar, double[] resArr) {
        for (int i = 0; i < inArr.length; i++) {
            resArr[i] = inArr[i] * scalar;
        }

        return resArr;
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

    public static void setMatrixRow (double[] aMat, double[] rowValuesToSet, int rowIdx, int dim) {
        System.arraycopy(rowValuesToSet, 0, aMat, rowIdx * dim, dim);
    }

    public static void getMatrixRow (double[] aMat, final int rowIdx, final int dim, final double[] rowValues) {
        for (int j = 0; j < dim; j++) {
            rowValues[j] = getMatrixEntry(aMat, rowIdx, j, dim);
        }
    }

    public static double getMatrixEntry(double[] aMat, final int rowIdx, final int colIdx, final int dim) {
        // row index i
        // column index j
        // dimension of the matrix (2D double array)
        // index in aMat[] = (i - 1) * dim + j
        return aMat[rowIdx * dim + colIdx];
    }

    /*
     * This method adds inVec to vecToAdd, i.e.
     * resVec <- inVec + vecToAdd
     */
    public static double[] vectorAdd (double[] inVec, double[] vecToAdd, double[] resVec) {
        for (int i = 0; i < inVec.length; i++) {
            resVec[i] = inVec[i] + vecToAdd[i];
        }
        return resVec;
    }

    /*
     * This method sets the value in a matrix
     *
     * row index i
     * column index j
     * dimension of the matrix (2D double array)
     * index in m[] = (i - 1) * dim + j
     */
    public static void setMatrixEntry(double[] m, final int i, final int j, final double value, final int dim){
        m[i * dim + j] = value;
    }

    /*
     * This method calculates vector dotProduct and matrix preMultiply, i.e.
     * scalar * inVec.transpose * inVec,
     * and returns a double value.
     */
    public static double vecTransScalarMultiply(double[] inVec, double scalar, int dim) {
        double res = 0.0;
        for (int i = 0; i < dim; i++) {
            res = res + inVec[i] * inVec[i];
        }
        return res * scalar;
    }

    /*
     * This method calculates the dot product of two vectors
     */
    public static double vectorDotMultiply(double[] aVec, double[] bVec) {
        double res = 0.0;

        for (int i = 0; i < aVec.length; i++) {
            res += aVec[i] * bVec[i];
        }
        return res;
    }

    /*
     * This method populates the matrix transpose
     * Requires a square matrix as inArr
     */
    public static double[] matrixTranspose (double [] inArr, int dim, double[] resArr) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                setMatrixEntry(resArr, j, i, getMatrixEntry(inArr, i, j, dim), dim);
            }
        }
        return resArr;
    }

    public static double[] matrixTranspose (double [] inArr, int nRow, int nCol, double[] resArr) {
        for (int i = 0; i < nRow; i++) {
            for (int j = 0; j < nCol; j++) {
                setMatrixEntry(resArr, j, i, getMatrixEntry(inArr, i, j, nCol), nRow);
            }
        }
        return resArr;
    }

    /*
     * This methods calculates product of inMat and matToMultiply.
     * The code that calls this function should already have verified
     * that (nColInMat == nRowMatToMultiply) is true.
     *
     * Also, note that nColInMat = nColResMat in matrix multiplication.
     */
    public static double[] matrixMultiply (double[] inMat, double[] matToMultiply, int nRowInMat, int nColInMat, double[] resMat) {
        for (int i = 0; i < nRowInMat; i ++) {
            for (int j = 0; j < nColInMat; j ++) {
                double sumTemp = 0.0;

                for (int k = 0; k < nColInMat; k++) {
                    sumTemp += getMatrixEntry(inMat, i, k, nColInMat) * getMatrixEntry(matToMultiply, k, j, nColInMat);
                }

                setMatrixEntry(resMat, i, j, sumTemp, nColInMat);
            }
        }

        return resMat;
    }
    /*
     * This methods calculates product of inMat and matToMultiply.
     * inMat has nRow and dim column number
     * matToMultiply has dim row number and nCol
     * resMat has nRow * nCol
     */
    public static void matrixMultiply (double[] inMat, double[] matToMultiply, int nRowInMat, int dim, int nColToMat, double[] resMat) {
        for (int i = 0; i < nRowInMat; i ++) {
            for (int j = 0; j < nColToMat; j ++) {
                double sumTemp = 0.0;

                for (int k = 0; k < dim; k++) {
                    sumTemp += MatrixUtilsContra.getMatrixEntry(inMat, i, k, dim) * MatrixUtilsContra.getMatrixEntry(matToMultiply, k, j, nColToMat);
                }

                MatrixUtilsContra.setMatrixEntry(resMat, i, j, sumTemp, nColToMat);
            }
        }
    }

    /*
     * This method multiplies matrix inRM by a double scalar.
     */
    public static double [] matrixScalarMultiply(double[] inArr, double scalar, int nColInMat, double[] resArr) {

        for (int i = 0; i < nColInMat; i++) {
            for (int j = 0; j < nColInMat; j++) {
                setMatrixEntry(resArr, i, j, getMatrixEntry(inArr, i, j, nColInMat) * scalar, nColInMat);
            }
        }
        return resArr;
    }

    /*
     * This method calculates
     * (1) product of inMat and matToMultiply.
     * (2) multiplication of the product and the scalar
     * The code that calls this function should already have verified
     * that (nColInMat == nRowMatToMultiply) is true.
     *
     * Also, note that nColInMat = nColResMat in matrix multiplication.
     */
    public static double[] matrixMultiplyScalar (double[] inMat, double[] matToMultiply, double scalar, int nRowInMat, int nColInMat, double[] resMat) {
        for (int i = 0; i < nRowInMat; i ++) {
            for (int j = 0; j < nColInMat; j ++) {
                double sumTemp = 0.0;

                for (int k = 0; k < nColInMat; k++) {
                    sumTemp += getMatrixEntry(inMat, i, k, nColInMat) * getMatrixEntry(matToMultiply, k, j, nColInMat);
                }

                setMatrixEntry(resMat, i, j, sumTemp * scalar, nColInMat);
            }
        }

        return resMat;
    }

    /*
     * This method calculates
     * (1) a row vector-matrix product,
     * (2) map multiply of the row vector,
     * i.e.,
     * resVec <- scalar * inVec * inRM
     *
     * Populates and returns resArr.
     */
    public static double[] matrixPreMapMultiply(double[] inArr, double[] inMatArr, double scalar, int nCol, int nRow, double[] resArr) {
        for (int i = 0; i < nCol; i++) {
            double sum = 0.0;

            for (int j = 0; j < nRow; j++) {
                sum = sum + inArr[j] * inMatArr[j * nCol + i];
            }

            resArr[i] = sum * scalar;
        }

        return resArr;
    }

    /*
     * This method calculates
     * (1) transpose of inMat
     * (1) a row vector-matrix product,
     * (2) map multiply of the row vector,
     * i.e.,
     * resVec <- scalar * inVec * inRM.transpose
     *
     * Populates and returns resArr.
     */
    public static double[] matrixTransPreMapMultiply(double[] inArr, double[] inMatArr, double scalar, int nCol, int nRow, double[] resArr) {
        for (int i = 0; i < nCol; i++) {
            double sum = 0.0;

            for (int j = 0; j < nRow; j++) {
                sum = sum + inArr[j] * inMatArr[i * nCol + j];
            }

            resArr[i] = sum * scalar;
        }

        return resArr;
    }

    /*
     * This methods calculates
     * (1) transpose of inMat
     * (2) product of inMat.transpose and matToMultiply.
     * i.e. inMat.transpose %*% matToMultiply
     * The code that calls this function should already have verified
     * that (nColInMat == nRowMatToMultiply) is true.
     *
     * Also, note that nColInMat = nColResMat in matrix multiplication.
     */
    public static double[] matrixTransMultiply (double[] inMat, double[] matToMultiply, int nRowInMat, int nColInMat, double[] resMat) {
        for (int i = 0; i < nRowInMat; i ++) {
            for (int j = 0; j < nColInMat; j ++) {
                double sumTemp = 0.0;

                for (int k = 0; k < nColInMat; k++) {
                    sumTemp += getMatrixEntry(inMat, k, i, nColInMat) * getMatrixEntry(matToMultiply, k, j, nColInMat);
                }

                setMatrixEntry(resMat, i, j, sumTemp, nColInMat);
            }
        }

        return resMat;
    }

    /*
     * This method calculates
     * scalar * (inMat + inMat.transpose)
     */
    public static double[] matrixTransAddScalar(double[] inMat, double scalar, int dim, double[] resMat) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                setMatrixEntry(resMat, i, j, (getMatrixEntry(inMat, i, j, dim) + getMatrixEntry(inMat, j, i, dim)) * scalar, dim);
            }
        }
        return resMat;
    }

    /*
     * This method calculates a row vector-matrix product, i.e.,
     * resVec <- inVec * inRM + vecToAdd
     *
     * Populates and returns resArr.
     */
    public static double[] matrixPreMultiplyAddVector (double[] inVec, double[] inMatArr, double[] vecToAdd, int nCol, int nRow, double[] resArr) {
        for (int i = 0; i < nCol; i++) {
            double sum = 0.0;

            for (int j = 0; j < nRow; j++) {
                sum = sum + inVec[j] * inMatArr[j * nCol + i];
            }

            resArr[i] = sum + vecToAdd[i];
        }

        return resArr;
    }

    /*
     * This method calculates
     * (1) transpose of inMat
     * (2) scalar product
     * i.e.
     * scalar * inMat.transpose
     */
    public static double[] matrixTransScalar(double[] inMat, double scalar, int dim, double[] resMat) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                setMatrixEntry(resMat, i, j, getMatrixEntry(inMat, j, i, dim) * scalar, dim);
            }
        }
        return resMat;
    }

    /*
     * This method subtract matrix rmToSub from inMat.
     *
     */
    public static double [] matrixSubtract(double[] inArr, double[] arrToSubtract, int dim, double[] resArr) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                setMatrixEntry(resArr, i, j, getMatrixEntry(inArr, i, j, dim) - getMatrixEntry(arrToSubtract, i, j, dim), dim);
            }
        }

        return resArr;
    }

    /*
     * This method subtracts vecToSub from inVec, i.e.
     * resVec <- inVec - vecToAdd
     */
    public static double[] vectorSubtract (double[] inVec, double[] vecToSub, double[] resVec) {
        for (int i = 0; i < inVec.length; i++) {
            resVec[i] = inVec[i] - vecToSub[i];
        }
        return resVec;
    }

    public static void populateMatrixArray(RealMatrix matrix, int nCol, int nRow, double[] array){
        for (int i = 0; i < nRow; i ++){
            System.arraycopy(matrix.getRow(i), 0, array, i * nCol, nCol);
        }
    }

    /*
     * This method calculates the multiplication of three matrices, i.e.
     * resMat = mat1 %*% mat2 %*% mat3
     *
     * NOTE: all matrices have to square with dimension "dim".
     */
    public static void matricesProduct(double[] mat1, double[] mat2, double[] mat3, int dim, double[] resMat) {
        for (int i = 0; i < dim; i++) {
            for (int l = 0; l < dim; l++) {
                double sum2 = 0.0;
                for (int k = 0; k < dim; k++) {
                    double sum1 = 0;
                    for (int j = 0; j < dim; j++) {
                        sum1 += MatrixUtilsContra.getMatrixEntry(mat1, i, j, dim) * MatrixUtilsContra.getMatrixEntry(mat2, j, k, dim);
                    }
                    sum2 += sum1 * MatrixUtilsContra.getMatrixEntry(mat3, k, l, dim);
                }
                MatrixUtilsContra.setMatrixEntry(resMat, i, l, sum2, dim);
            }
        }
    }

    /*
     * This method calculates the multiplication of three matrices, i.e.
     * resMat = mat1 %*% mat2 %*% t(mat1) + mat3
     *
     * NOTE: all matrices have to square with dimension "dim".
     */
    public static void matricesTransProductAdd(double[] mat1, double[] mat2, double[] mat3, int dim, double[] resMat){
        for (int i = 0; i < dim; i++) {
            for (int l = 0; l < dim; l++) {
                double sum2 = 0.0;
                for (int k = 0; k < dim; k++) {
                    double sum1 = 0;
                    for (int j = 0; j < dim; j++) {
                        sum1 += MatrixUtilsContra.getMatrixEntry(mat1, i, j, dim) * MatrixUtilsContra.getMatrixEntry(mat2, j, k, dim);
                    }
                    sum2 += sum1 * MatrixUtilsContra.getMatrixEntry(mat1, l, k, dim);
                }
                MatrixUtilsContra.setMatrixEntry(resMat, i, l, sum2 + MatrixUtilsContra.getMatrixEntry(mat3, i,l,dim), dim);
            }
        }

    }

}
