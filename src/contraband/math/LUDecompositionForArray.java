package contraband.math;

/*
 *  This class rewrites "org.apache.commons LUDecomposition",
 *  by using 1D double array.
 */

import org.apache.commons.math3.util.FastMath;

public class LUDecompositionForArray {
    private static final double DEFAULT_TOO_SMALL = 1e-11;

    public static void ArrayLUDecomposition(double[] matrix, double[] lu, int[] pivot, boolean[] evensingular, final int m) {

        System.arraycopy(matrix,0, lu, 0, m * m);

        // Initialize permutation array and parity
        for (int row = 0; row < m; row++) {
            pivot[row] = row;
        }


        evensingular[0] = true;
        evensingular[1] = false;

        // Loop over columns
        for (int col = 0; col < m; col++) {

            // upper
            for (int row = 0; row < col; row++) {
                double sum = MatrixUtilsContra.getMatrixEntry(lu, row, col, m);
                for (int i = 0; i < row; i++) {
                    sum -= MatrixUtilsContra.getMatrixEntry(lu, row, i, m) * MatrixUtilsContra.getMatrixEntry(lu, i, col, m);
                }
                MatrixUtilsContra.setMatrixEntry(lu, row, col, sum, m);
            }

            // lower
            int max = col; // permutation row
            double largest = Double.NEGATIVE_INFINITY;
            for (int row = col; row < m; row++) {
                double sum = MatrixUtilsContra.getMatrixEntry(lu, row, col, m);
                for (int i = 0; i < col; i++) {
                    sum -= MatrixUtilsContra.getMatrixEntry(lu, row, i, m) * MatrixUtilsContra.getMatrixEntry(lu, i, col, m);
                }
                //luRow[col] = sum;
                MatrixUtilsContra.setMatrixEntry(lu, row, col, sum, m);

                // maintain best permutation choice
                if (FastMath.abs(sum) > largest) {
                    largest = FastMath.abs(sum);
                    max = row;
                }
            }

            // Singularity check
            //if (FastMath.abs(lu[max][col]) < DEFAULT_TOO_SMALL) {
            if (FastMath.abs(MatrixUtilsContra.getMatrixEntry(lu, max, col, m)) < DEFAULT_TOO_SMALL) {
                //singular = true;
                //return new boolean[] {even, singular};
                evensingular[1] = true;
            }

            // Pivot if necessary
            if (max != col) {
                double tmp = 0;
                for (int i = 0; i < m; i++) {
                    tmp = MatrixUtilsContra.getMatrixEntry(lu, max, i, m);
                    MatrixUtilsContra.setMatrixEntry(lu, max, i, MatrixUtilsContra.getMatrixEntry(lu, col, i,m), m);
                    MatrixUtilsContra.setMatrixEntry(lu, col, i, tmp, m);
                }


                int temp = pivot[max];
                pivot[max] = pivot[col];
                pivot[col] = temp;
                //even = !even;
                evensingular[0] = !evensingular[0];
            }

            // Divide the lower elements by the "winning" diagonal elt.
            //final double luDiag = lu[col][col];
            final double luDiag = MatrixUtilsContra.getMatrixEntry(lu, col, col, m);
            for (int row = col + 1; row < m; row++) {
                //lu[row][col] /= luDiag;
                MatrixUtilsContra.setMatrixEntry(lu, row, col, MatrixUtilsContra.getMatrixEntry(lu, row, col, m) / luDiag, m);
            }
        }

    }

    public static double getDeterminant(double[] lu, int m, boolean[] evensingular) {
        boolean even = evensingular[0];
        boolean singular = evensingular[1];

        if (singular) {
            return 0;
        } else {
            double determinant = even ? 1 : -1;
            for (int i = 0; i < m; i++) {
                determinant = determinant * MatrixUtilsContra.getMatrixEntry(lu, i, i, m);
            }
            return determinant;
        }
    }

    /**
     * Get the inverse of the decomposed matrix.
     */
    public static void populateInverseMatrix(double[] lu, int[] pivot, double [] identity, boolean singular, int m, double[] inv) {

        if (singular) {
            throw new RuntimeException("Singular matrix!");
        }

        //final double[] bpRow = new double [m];

        for (int row = 0; row < m; row++) {
            //final double[] bpRow = inv[row];
            final int pRow = pivot[row];
            for (int col = 0; col < m; col++) {
                //bpRow[col] = identity[pRow][col];
                MatrixUtilsContra.setMatrixEntry(inv, row, col, MatrixUtilsContra.getMatrixEntry(identity, pRow, col, m), m);
            }

        }

        // Solve LY = b
        for (int col = 0; col < m; col++) {
            //final double[] bpCol = inv[col];

            for (int i = col + 1; i < m; i++) {
                //final double[] bpI = inv[i];

                //final double luICol = lu[i][col];
                final double luICol = MatrixUtilsContra.getMatrixEntry(lu, i, col, m);

                for (int j = 0; j < m; j++) {
                    //bpI[j] -= bpCol[j] * luICol;

                    // bpCol[j] = inv[j][col]
                    // bpI[j] = inv[i][j]
                    double bpIj = MatrixUtilsContra.getMatrixEntry(inv, i, j, m);
                    double bpColj = MatrixUtilsContra.getMatrixEntry(inv, col, j, m);
                    MatrixUtilsContra.setMatrixEntry(inv, i, j, bpIj - bpColj * luICol, m);
                }
            }
        }

        // Solve UX = Y
        for (int col = m - 1; col >= 0; col--) {
            //final double[] bpCol = inv[col];

            //final double luDiag = lu[col][col];
            final double luDiag = MatrixUtilsContra.getMatrixEntry(lu, col, col, m);

            for (int j = 0; j < m; j++) {
                //bpCol[j] /= luDiag;
                double boColj = MatrixUtilsContra.getMatrixEntry(inv, col, j, m);
                MatrixUtilsContra.setMatrixEntry(inv, col, j, boColj / luDiag, m);


            }
            for (int i = 0; i < col; i++) {
                //final double[] bpI = inv[i];

                //final double luICol = lu[i][col];
                final double luICol = MatrixUtilsContra.getMatrixEntry(lu, i, col, m);

                for (int j = 0; j < m; j++) {
                    //bpI[j] -= bpCol[j] * luICol;
                    double bpIj = MatrixUtilsContra.getMatrixEntry(inv, i, j, m);
                    double bpColj = MatrixUtilsContra.getMatrixEntry(inv, col, j, m);
                    MatrixUtilsContra.setMatrixEntry(inv, i, j, bpIj-bpColj*luICol,m);


                }
            }

        }

    }
}
