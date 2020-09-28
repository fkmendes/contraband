package contraband.utils;

import contraband.math.MatrixUtils;

import static java.lang.Math.exp;

// author: Pau Bravo
// title: mvMORPH: An R package for fitting multivariate evolutionary models to morphometric data 
// Equation 9
public class EBUtils {
	
	public static void computeEBtMat(double[][] gMat, double[][] Sigma, double[][] ebMat) {
		MatrixUtils.kronecker(gMat, Sigma, ebMat);
	}
	
	// Computation of the g speed matrix for EB model ( (exp(g * tMat) - 1)/g )
    public static void computeGSpeedMat(double gSpeed, double[][] phyloTMat, double[][] gMat) { // gMat must have dimensions n x n like phyloTMat
        int n = phyloTMat.length;

        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                gMat[row][col] = (exp(gSpeed * phyloTMat[row][col]) - 1) / gSpeed;
            }
        }
    }
		
}
