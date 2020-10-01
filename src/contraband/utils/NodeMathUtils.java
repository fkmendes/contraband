package contraband.utils;

import contraband.math.MatrixUtilsContra;
import org.apache.commons.math3.util.FastMath;

public class NodeMathUtils {
    /*
     * The following two methods are used to populate the Sigma matrix by using correlations and one rate for all traits.
     * Note: "traitRateMatrix" <- UpperSigma * t(UpperSigma)
     */
    public static void populateUpperTriangularMat(double variance, double[] covariance, int nTraits, double[] upperTriangularMat) {
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(upperTriangularMat, i, i, variance, nTraits);
            for (int j = i + 1; j < nTraits; j++) {
                double cov = variance * covariance[k];
                MatrixUtilsContra.setMatrixEntry(upperTriangularMat,i, j, cov, nTraits);
                k++;
            }
        }
    }

    public static void populateTraitRateMatrix(double variance, double[] covariance, double[] sigma, double [] transSigma, int nTraits, double[] traitRateMatrix) {
        populateUpperTriangularMat(variance, covariance, nTraits, sigma);

        MatrixUtilsContra.vectorTranspose(sigma, nTraits, transSigma);

        MatrixUtilsContra.matrixMultiply(sigma, transSigma, nTraits, traitRateMatrix);
    }

    /*
     * The following two methods are used to populate the Sigma matrix by using correlations
     * and each trait has its own rate.
     * Note: "traitRateMatrix" <- UpperSigma * t(UpperSigma)
     */
    public static void populateUpperTriangularMat(double[] variance, double[] covariance, int nTraits, double[] upperTriangularMat) {
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(upperTriangularMat, i, i, variance[i], nTraits);
            for (int j = i + 1; j < nTraits; j++) {
                double cov = FastMath.sqrt(variance[i]) * FastMath.sqrt(variance[j]) * covariance[k];
                MatrixUtilsContra.setMatrixEntry(upperTriangularMat,i, j, cov, nTraits);
                k++;
            }
        }
    }

    /*
     * This method populate UpperTriangular matrix for the traitRateMatrix
     * by using the input variance and covariance.
     *
     */
    public static void populateCoUpperTriangularMat(double[] variance, double[] covariance, int nTraits, double[] upperTriangularMat) {
        int k = 0;
        int m = 0;
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(upperTriangularMat, i, i, variance[m], nTraits);
            m++;
            for (int j = i + 1; j < nTraits; j++) {
                MatrixUtilsContra.setMatrixEntry(upperTriangularMat,i, j, covariance[k], nTraits);
                k++;
            }
        }
    }

    public static void populateTraitRateMatrix(double[] variance, double[] covariance, double[] sigma, double [] transSigma, int nTraits, double[] traitRateMatrix, boolean coEstimate) {
        if(coEstimate) {
            populateCoUpperTriangularMat(variance, covariance, nTraits, sigma);
        } else {
            populateUpperTriangularMat(variance, covariance, nTraits, sigma);
        }

        MatrixUtilsContra.vectorTranspose(sigma, nTraits, transSigma);

        MatrixUtilsContra.matrixMultiply(sigma, transSigma, nTraits, traitRateMatrix);
    }

    public static void populateTraitRateMatrixDirectly(double[] sigmasq, double[] rho, int nTraits, double[] traitRateMatrix) {
        // diagonal = sigma^2
        // off-diagonal = sigma_i * sigma_j * rho_ij
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            double var = sigmasq[i];
            MatrixUtilsContra.setMatrixEntry(traitRateMatrix, i, i, var, nTraits);
            for (int j = i + 1; j < nTraits; j++) {
                double cov = FastMath.sqrt(sigmasq[i]) * FastMath.sqrt(sigmasq[j]) * rho[k];
                MatrixUtilsContra.setMatrixEntry(traitRateMatrix, i, j, cov, nTraits);
                MatrixUtilsContra.setMatrixEntry(traitRateMatrix, j, i, cov, nTraits);
                k++;
            }
        }
    }

    public static void populateTraitRateMatrixDirectly(double sigmasq, double[] rho, int nTraits, double[] traitRateMatrix) {
        // R = sigma * Rho
        int k = 0;
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(traitRateMatrix, i, i, sigmasq, nTraits);
            for (int j = i + 1; j < nTraits; j++) {
                double cov = sigmasq * rho[k];
                MatrixUtilsContra.setMatrixEntry(traitRateMatrix,i, j, cov, nTraits);
                MatrixUtilsContra.setMatrixEntry(traitRateMatrix,j, i, cov, nTraits);
                k++;
            }
        }
    }
}
