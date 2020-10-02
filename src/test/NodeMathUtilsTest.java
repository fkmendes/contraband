package test;

import contraband.utils.NodeMathUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.junit.Assert;
import org.junit.Test;

public class NodeMathUtilsTest {

    /*
     * In R
     *
     * a = matrix(c(1.0, 0.3, 0.2,
     *              0.0, 1.0, -0.1,
     *              0.0, 0.0, 1.0), ncol = 3, byrow = T)
     *
     * res = 2.5 * a
     */
    @Test
    public void testPopulateUpperTriangularMatForOneRate () {
        int nTraits = 3;
        double sigmasq = 2.5;
        double[] correlations = new double[] {0.3, 0.2, -0.1};

        double[] resMat = new double[nTraits * nTraits];

        NodeMathUtils.populateUpperTriangularMat(sigmasq, correlations, nTraits, resMat);

        Assert.assertArrayEquals(new Double[] {2.5, 0.75, 0.5, 0.0, 2.5, -0.25, 0.0, 0.0, 2.5}, ArrayUtils.toObject(resMat));
    }

    /*
     * a = matrix(c(1.0, 0.3, 0.2,
     *              0.0, 1.0, -0.1,
     *              0.0, 0.0, 1.0), ncol = 3, byrow = T)
     *
     * diag(a) = c(25, 9, 16)
     * a[1,2] = sqrt(25) * sqrt(9) * 0.3
     * a[1,3] = sqrt(25) * sqrt(16) * 0.3
     * a[2,3] = sqrt(9) * sqrt(16) * (-0.1)
     *
     */
    @Test
    public void testPopulateUpperTriangularMatForMultipleRates () {
        int nTraits = 3;
        double[] sigmasq = new double [] {25, 9, 16};
        double[] correlations = new double[] {0.3, 0.2, -0.1};

        double[] resMat = new double[nTraits * nTraits];

        NodeMathUtils.populateUpperTriangularMat(sigmasq, correlations, nTraits, resMat);

        Assert.assertArrayEquals(new Double[] {25.0, 4.5, 4.0, 0.0, 9.0, -1.2000000000000002, 0.0, 0.0, 16.0}, ArrayUtils.toObject(resMat));
    }

    /*
     * In R
     *
     * a = matrix(c(1.0, 0.3, 0.2,
     *              0.3, 1.0, -0.1,
     *              0.2, -0.1, 1.0), ncol = 3, byrow = T)
     *
     * res = 2.5 * a
     */
    @Test
    public void testPopulateTraitRateMatrixDirectlyForOneRate () {
        int nTraits = 3;
        double sigmasq = 2.5;
        double[] correlations = new double[] {0.3, 0.2, -0.1};

        double[] resMat = new double[nTraits * nTraits];

        NodeMathUtils.populateTraitRateMatrixDirectly(sigmasq, correlations, nTraits, resMat);

        Assert.assertArrayEquals(new Double[] {2.5, 0.75, 0.5, 0.75, 2.5, -0.25, 0.5, -0.25, 2.5}, ArrayUtils.toObject(resMat));
    }

    /*
     * In R
     *
     * a = matrix(c(1.0, 0.3, 0.2,
     *              0.3, 1.0, -0.1,
     *              0.2, -0.1, 1.0), ncol = 3, byrow = T)
     *
     * diag(a) = c(25, 9, 16)
     * a[1,2] = sqrt(25) * sqrt(9) * 0.3
     * a[2,1] = a[1,2]
     * a[1,3] = sqrt(25) * sqrt(16) * 0.2
     * a[3,1] = a[1,3]
     * a[2,3] = sqrt(9) * sqrt(16) * (-0.1)
     * a[3,2] = a[2,3]
     */
    @Test
    public void testPopulateTraitRateMatrixDirectlyForMultipleRates () {
        int nTraits = 3;
        double[] sigmasq = new double [] {25, 9, 16};
        double[] correlations = new double[] {0.3, 0.2, -0.1};

        double[] resMat = new double[nTraits * nTraits];

        NodeMathUtils.populateTraitRateMatrixDirectly(sigmasq, correlations, nTraits, resMat);

        Assert.assertArrayEquals(new Double[] {25.0, 4.5, 4.0, 4.5, 9.0, -1.2000000000000002, 4.0, -1.2000000000000002, 16.0}, ArrayUtils.toObject(resMat));
    }

    /*
     * In R
     *
     * a = matrix(c(1.0, 0.3, 0.2,
     *              0.0, 1.0, -0.1,
     *              0.0, 0.0, 1.0), ncol = 3, byrow = T)
     *
     * b = 2.5 * a
     *
     * res = b %*% t(b)
     */
    @Test
    public void testPopulateTraitRateMatrixForOneRate () {
        int nTraits = 3;
        double sigmasq = 2.5;
        double[] correlations = new double[] {0.3, 0.2, -0.1};

        double[] resMat = new double[nTraits * nTraits];

        double[] upperMatrix = new double[nTraits * nTraits];
        double[] transMatrix = new double[nTraits * nTraits];

        NodeMathUtils.populateTraitRateMatrix(sigmasq, correlations, upperMatrix, transMatrix, nTraits, resMat);

        Assert.assertArrayEquals(new Double[] {7.0625, 1.75, 1.25, 1.75, 6.3125, -0.625, 1.25, -0.625, 6.25}, ArrayUtils.toObject(resMat));
    }

    /*
     * In R
     *
     * res = matrix(c(25.0, 0.3, 0.2,
     *                0.0, 9.0, -0.1,
     *                0.0, 0.0, 16.0), ncol = 3, byrow = T)
     */
    @Test
    public void testPopulateCoUpperTriangularMat () {
        int nTraits = 3;
        double[] sigmasq = new double [] {25, 9, 16};
        double[] correlations = new double[] {0.3, 0.2, -0.1};

        double[] resMat = new double[nTraits * nTraits];

        NodeMathUtils.populateCoUpperTriangularMat(sigmasq, correlations, nTraits, resMat);

        Assert.assertArrayEquals(new Double[] {25.0, 0.3, 0.2, 0.0, 9.0, -0.1, 0.0, 0.0, 16.0}, ArrayUtils.toObject(resMat));
    }

    /*
     * In R
     *
     * a = matrix(c(25.0, 0.3, 0.2,
     *                0.0, 9.0, -0.1,
     *                0.0, 0.0, 16.0), ncol = 3, byrow = T)
     * res1 = a %*% t(a)
     *
     * a[1,2] = sqrt(25) * sqrt(9) * 0.3
     * a[1,3] = sqrt(25) * sqrt(16) * 0.2
     * a[2,3] = sqrt(9) * sqrt(16) * (-0.1)
     *
     * res2 = a %*% t(a)
     */
    @Test
    public void testPopulateTraitRateMatrixForMultipleRates () {
        int nTraits = 3;
        double[] sigmasq = new double [] {25, 9, 16};
        double[] correlations = new double[] {0.3, 0.2, -0.1};

        double[] resMat1 = new double[nTraits * nTraits];
        double[] resMat2 = new double[nTraits * nTraits];

        double[] upperMatrix = new double[nTraits * nTraits];
        double[] transMatrix = new double[nTraits * nTraits];

        NodeMathUtils.populateTraitRateMatrix(sigmasq, correlations, upperMatrix, transMatrix, nTraits, resMat1, true);
        Assert.assertArrayEquals(new Double[] {625.13, 2.6799999999999997, 3.2, 2.6799999999999997, 81.01, -1.6, 3.2, -1.6, 256.0}, ArrayUtils.toObject(resMat1));

        NodeMathUtils.populateTraitRateMatrix(sigmasq, correlations, upperMatrix, transMatrix, nTraits, resMat2, false);
        Assert.assertArrayEquals(new Double[] {661.25, 35.7, 64.0, 35.7, 82.44, -19.200000000000003, 64.0, -19.200000000000003, 256.0}, ArrayUtils.toObject(resMat2));
    }
}
