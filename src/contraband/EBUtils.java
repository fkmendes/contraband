package contraband;

import static java.lang.Math.exp;

// author: Clavel, Julien and Escarguel, Gilles and Merceron, Gildas
// title: mvMORPH: An R package for fitting multivariate evolutionary models to morphometric data 
//Equation 9
public class EBUtils {
	
	
	public static void computeEBtMat(double[][] gMat, double[][] Sigma, double[][] EBMat) {
		
		GeneralUtils.kronecker(gMat, Sigma, EBMat);
	}
	
	// Computation of the g speed matrix for EB model ( (exp(g * tMat) - 1)/g )
    public static void Gchunk(double gspeed, double[][] phyloMat, double[][] gMat) { // gMat must have dimensions n x n like phyloMat

        int n = phyloMat.length;

        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                gMat[row][col] = (exp(gspeed * phyloMat[row][col]) - 1) / gspeed;
            }
        }
    }
		
}
