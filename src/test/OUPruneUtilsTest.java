package test;

import contraband.math.MatrixUtilsContra;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;

public class OUPruneUtilsTest {
    final static double EPSILON = 1e-4;

    /*
	 * In R:
	 * a <- matrix(c(17.0, 21.0, 25.0, 29.0,
	 *               18.0, 22.0, 26.0, 30.0,
	 *               19.0, 23.0, 27.0, 31.0,
	 *               20.0, 24.0, 28.0, 32.0), ncol=4, byrow=T)
	 *
	 * res = 0.25 * a
	 *
     *     [,1] [,2] [,3] [,4]
     *[1,] 4.25 5.25 6.25 7.25
     *[2,] 4.50 5.50 6.50 7.50
     *[3,] 4.75 5.75 6.75 7.75
     *[4,] 5.00 6.00 7.00 8.00
     */
    @Test
    public void testMatrixScalarMultiply() {
        int ncol = 4;
        double[] aMat = new double[]{
                17.0, 21.0, 25.0, 29.0,
                18.0, 22.0, 26.0, 30.0,
                19.0, 23.0, 27.0, 31.0,
                20.0, 24.0, 28.0, 32.0
        };
        double[] resMat = new double[ncol * ncol];

        double scalar = 0.25;

        Double[] resMatDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixScalarMultiply(aMat, scalar, ncol, resMat));

        Assert.assertArrayEquals(new Double[]{4.25, 5.25, 6.25, 7.25, 4.5, 5.5, 6.5, 7.5, 4.75, 5.75, 6.75, 7.75, 5.0, 6.0, 7.0, 8.0}, resMatDouble);
    }

    /*
     * In R:
     * a <- matrix(c(1.0, 5.0, 9.0, 13.0,
     *               2.0, 6.0, 10.0, 14.0,
     *               3.0, 7.0, 11.0, 15.0,
     *               4.0, 8.0, 12.0, 16.0), ncol=4, byrow=T)
     * b <- matrix(c(17.0, 21.0, 25.0, 29.0,
     *               18.0, 22.0, 26.0, 30.0,
     *               19.0, 23.0, 27.0, 31.0,
     *               20.0, 24.0, 28.0, 32.0), ncol=4, byrow=T)
     *
     * 0.5 * a %*% b
     *       [,1] [,2] [,3] [,4]
     * [1,]  269  325  381  437
     * [2,]  306  370  434  498
     * [3,]  343  415  487  559
     * [4,]  380  460  540  620
     */
    @Test
    public void testMatrixMultiplyScalar() {
        int nRowInMat = 4;
        int nColInMat = 4;
        int length = nRowInMat * nColInMat;

        double[] aMat = new double[] {
                1.0, 5.0, 9.0, 13.0,
                2.0, 6.0, 10.0, 14.0,
                3.0, 7.0, 11.0, 15.0,
                4.0, 8.0, 12.0, 16.0
        };

        double[] bMat = new double[] {
                17.0, 21.0, 25.0, 29.0,
                18.0, 22.0, 26.0, 30.0,
                19.0, 23.0, 27.0, 31.0,
                20.0, 24.0, 28.0, 32.0
        };

        double[] resMat = new double[length];

        double scalar = 0.5;
        Double[] resMatDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixMultiplyScalar(aMat, bMat, scalar, nRowInMat, nColInMat, resMat));

        Assert.assertArrayEquals(new Double[] {
                269.0, 325.0, 381.0, 437.0, 306.0, 370.0, 434.0, 498.0, 343.0, 415.0, 487.0, 559.0, 380.0, 460.0, 540.0, 620.0 }, resMatDouble);
    }

    /*
     * In R:
     *
     * a <- matrix(c(1,2,3,4,5,6,7,8,9), nrow=3)
     * b <- c(1,2,3)
     * t(b)
     *       [,1] [,2] [,3]
     * [1,]    1    2    3
     *
     * 0.5 * t(b) %*% a
     *       [,1] [,2] [,3]
     * [1,]   7   16   25
     */
    @Test
    public void testMatrixPreMapMultiply() {
        double[] inArr = new double[] { 1, 2, 3 };
        double[] inMatArr = new double[] { 1, 4, 7, 2, 5, 8, 3, 6, 9 };
        double[] resArr = new double[3];

        int nCol = 3;
        int nRow = 3;

        double scalar = 0.5;

        Double[] resArrDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixPreMapMultiply(inArr, inMatArr, scalar, nCol, nRow, resArr));

        Assert.assertArrayEquals(resArrDouble, new Double[] { 7.0, 16.0, 25.0 });
    }

    /*
     * In R:
     *
     * a <- matrix(1:16, nrow=4)
     * b <- c(1,2,3,4)
     * t(b)
     *       [,1] [,2] [,3] [,4]
     * [1,]    1    2    3   4
     *
     * 0.5 * t(b) %*% t(a)
     *       [,1] [,2] [,3]  [,4]
     * [1,]   45   50   55   60
     */
    @Test
    public void testMatrixTransPreMapMultiply() {
        double[] inArr = new double[] { 1, 2, 3, 4 };
        double[] inMatArr = new double[] { 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16 };
        double[] resArr = new double[4];

        int nCol = 4;
        int nRow = 4;

        double scalar = 0.5;

        Double[] resArrDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixTransPreMapMultiply(inArr, inMatArr, scalar, nCol, nRow, resArr));

        Assert.assertArrayEquals(resArrDouble, new Double[] { 45.0, 50.0, 55.0, 60.0 });
    }

    /*
     * In R:
     *
     * a <- matrix(1:16, nrow=4)
     * b <- matrix(0.2 * (1:16), nrow=4)
     *
     *  t(a) %*% b
     *         [,1] [,2] [,3]  [,4]
     *  [1,]    6 14.0  22.0  30.0
     *  [2,]   14 34.8  55.6  76.4
     *  [3,]   22 55.6  89.2 122.8
     *  [4,]   30 76.4 122.8 169.2
     */
    @Test
    public void testMatrixTransMultiply() {
        double[] inMatArr = new double[] { 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16 };
        double[] matToMul = new double[] { 0.2, 1, 1.8, 2.6, 0.4, 1.2, 2, 2.8, 0.6, 1.4, 2.2, 3, 0.8, 1.6, 2.4, 3.2 };
        double[] resMat = new double[16];

        int nCol = 4;
        int nRow = 4;

        Double[] resArrDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixTransMultiply(inMatArr, matToMul, nCol, nRow, resMat));

        Assert.assertArrayEquals(resArrDouble, new Double[] { 6.0, 14.0, 22.0, 30.0, 14.0, 34.8, 55.60000000000001, 76.4, 22.0, 55.6, 89.2, 122.80000000000001, 30.000000000000004, 76.4, 122.80000000000001, 169.2 });
    }

     /*
     * In R:
     *
     * a <- matrix(c(2.550531, 2.350460, 1.242473, 1.534543,
     *              1.702243, 2.455985, 2.423319, 1.267112,
     *              2.077579, 2.054758, 2.861430, 2.958108,
     *              2.274022, 1.945004, 1.524684, 1.960117), nrow=4, byrow=T)
     *
     * 3.2 * (a + t(a))
     *
     *       [,1]     [,2]     [,3]     [,4]
     * [1,] 16.32340 12.96865 10.62417 12.18741
     * [2,] 12.96865 15.71830 14.32985 10.27877
     * [3,] 10.62417 14.32985 18.31315 14.34493
     * [4,] 12.18741 10.27877 14.34493 12.54475
     *
     */
    @Test
    public void testMatrixScalarAdd() {
        double[] inMatArr = new double[] {2.550531, 2.350460, 1.242473, 1.534543,
                                          1.702243, 2.455985, 2.423319, 1.267112,
                                          2.077579, 2.054758, 2.861430, 2.958108,
                                          2.274022, 1.945004, 1.524684, 1.960117};
        double[] resMat = new double[16];
        double scalar = 3.2;

        int nCol = 4;

        Double[] resArrDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixTransAddScalar(inMatArr, scalar, nCol, resMat));

        Assert.assertArrayEquals(resArrDouble, new Double[] {16.3233984, 12.9686496, 10.6241664, 12.187408, 12.9686496, 15.718304000000002, 14.329846400000001, 10.278771200000001, 10.6241664, 14.329846400000001, 18.313152, 14.3449344, 12.187408, 10.278771200000001, 14.3449344, 12.5447488});
    }

    /*
     * In R:
     *
     * a <- matrix(1:16, nrow=4)
     * b <- c(1:4)
     * c <- c(2:5)
     *
     * c + b %*% a
     *       [,1] [,2] [,3] [,4]
     * [1,]   32   73  114  155
     *
     */
    @Test
    public void testMatrixPreMultiplyAddVector() {
        double[] inMatArr = new double[]{1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16};
        double[] inVec = new double[]{1, 2, 3, 4};
        double[] vecToAdd = new double[]{2, 3, 4, 5};

        double[] resMat = new double[4];
        int nCol = 4;
        int nRow = 4;

        Double[] resArrDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixPreMultiplyAddVector(inVec, inMatArr, vecToAdd, nCol, nRow, resMat));

        Assert.assertArrayEquals(resArrDouble, new Double[] {32.0, 73.0, 114.0, 155.0});
    }

    /*
     * In R:
     *
     * a <- matrix(c(2.550531, 2.350460, 1.242473, 1.534543,
     *              1.702243, 2.455985, 2.423319, 1.267112,
     *              2.077579, 2.054758, 2.861430, 2.958108,
     *              2.274022, 1.945004, 1.524684, 1.960117), nrow=4, byrow=T)
     *
     * 0.2 * t(a)
     *
     */
    @Test
    public void testMatrixTransScalar() {
        double[] inMatArr = new double[] {2.550531, 2.350460, 1.242473, 1.534543,
                1.702243, 2.455985, 2.423319, 1.267112,
                2.077579, 2.054758, 2.861430, 2.958108,
                2.274022, 1.945004, 1.524684, 1.960117};
        double[] resMat = new double[16];
        double scalar = 0.2;

        int nCol = 4;

        Double[] resArrDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixTransScalar(inMatArr, scalar, nCol,resMat));

        Assert.assertArrayEquals(resArrDouble, new Double[] {0.5101062, 0.3404486, 0.41551580000000005, 0.4548044, 0.470092, 0.49119700000000005, 0.41095160000000003, 0.38900080000000004, 0.2484946, 0.4846638, 0.572286, 0.3049368, 0.30690860000000003, 0.2534224, 0.5916216000000001, 0.3920234});
    }

    /*
     * In R:
     *
     * a <- matrix(1:16, nrow=4)
     * b <- matrix(0.5 * (5:20), nrow=4)
     *
     *  a - b
     *        [,1] [,2] [,3] [,4]
     *   [1,] -1.5  0.5  2.5  4.5
     *   [2,] -1.0  1.0  3.0  5.0
     *   [3,] -0.5  1.5  3.5  5.5
     *   [4,]  0.0  2.0  4.0  6.0
     *
     */
    @Test
    public void testMatrixSubtract(){
        double[] inMatArr = new double[]{1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16};
        double[] matToSub = new double[]{2.5, 4.5, 6.5, 8.5, 3, 5, 7, 9, 3.5, 5.5, 7.5, 9.5, 4, 6, 8, 10};

        double[] resMat = new double[16];

        int nCol = 4;

        Double[] resArrDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixSubtract(inMatArr, matToSub, nCol, resMat));

        Assert.assertArrayEquals(resArrDouble, new Double[] {-1.5, 0.5, 2.5, 4.5, -1.0, 1.0, 3.0, 5.0, -0.5, 1.5, 3.5, 5.5, 0.0, 2.0, 4.0, 6.0});
    }

    /*
     * In R:
     *
     * a <- c(1.5, 5.8, 9.6, 13.1)
     * b <- c(6.3, 2.3, 4.3, 1.2)
     *
     * a - b
     *      [,1] [,2] [,3] [,4]
     * [1,] -4.8  3.5  5.3 11.9
     */
    @Test
    public void testVectorSubtract() {
        double[] inVec = new double[]{1.5, 5.8, 9.6, 13.1};
        double[] vecToSub = new double[] {6.3, 2.3, 4.3, 1.2};

        double[] resVec = new double[4];

        Double[] resArrDouble = ArrayUtils.toObject(MatrixUtilsContra.vectorSubtract(inVec, vecToSub, resVec));
        Assert.assertArrayEquals(resArrDouble, new Double[] {-4.8, 3.5, 5.3, 11.9});
    }

    @Test
    public void testPopulateMatrixArray() {
        int nCol = 8;
        int nRow = 7;
        RealMatrix matRM = new Array2DRowRealMatrix(new double [][]
                {{5.555278, 3.637867, 3.826537, 5.580623, 4.233251, 4.265012, 5.668205, 5.960273},
                        {7.102671, 3.008732, 3.668539, 7.797644, 7.149419, 4.708736, 6.938855, 6.028421},
                        {8.244949, 2.804517, 5.590706, 4.507473, 3.285521, 5.566429, 1.413961, 6.178505},
                        {6.623658, 7.340270, 2.815117, 5.165536, 6.636394, 8.104003, 4.286587, 3.882813},
                        {5.732623, 1.996666, 6.853803, 5.431592, 1.740781, 7.074177, 6.466281, 2.540964},
                        {6.371464, 8.338275, 4.371867, 8.827669, 7.309662, 5.661179, 8.103298, 2.239197},
                        {1.422262, 1.983072, 4.297063, 4.046804, 5.934354, 7.656025, 5.650181, 2.041209}
                });

        double[] array = new double[nCol * nRow];

        MatrixUtilsContra.populateMatrixArray(matRM, nCol, nRow, array);
        Double[] resArrDouble = ArrayUtils.toObject(array);
        Assert.assertArrayEquals(resArrDouble, new Double[] {5.555278, 3.637867, 3.826537, 5.580623, 4.233251, 4.265012, 5.668205, 5.960273,
                7.102671, 3.008732, 3.668539, 7.797644, 7.149419, 4.708736, 6.938855, 6.028421,
                8.244949, 2.804517, 5.590706, 4.507473, 3.285521, 5.566429, 1.413961, 6.178505,
                6.623658, 7.340270, 2.815117, 5.165536, 6.636394, 8.104003, 4.286587, 3.882813,
                5.732623, 1.996666, 6.853803, 5.431592, 1.740781, 7.074177, 6.466281, 2.540964,
                6.371464, 8.338275, 4.371867, 8.827669, 7.309662, 5.661179, 8.103298, 2.239197,
                1.422262, 1.983072, 4.297063, 4.046804, 5.934354, 7.656025, 5.650181, 2.041209});
    }



}
