package contraband;

import org.apache.commons.math3.linear.RealMatrix;

public class OUUtils {
	
	/*
	 * Note: sigma^2 does not go in because it is handled by MVNUtils
	 */
	public static RealMatrix computeOUTMatOneTrait(int n, double alpha, RealMatrix tMat, RealMatrix ouTMat, boolean rootIsFixed) {
		
		double cellValue;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				cellValue = (tMat.getEntry(i, i) + tMat.getEntry(j, j) - 2*tMat.getEntry(i, j));
				
				if (rootIsFixed) {
					cellValue *= (1 - Math.exp(-2 * alpha * tMat.getEntry(i, j)));
				}
				
				ouTMat.setEntry(i, j, cellValue);  // exponent part
			}
		}
			
		ouTMat = tMat.scalarMultiply(1 / (2*alpha)); // divide by 2*alpha
		
		return ouTMat;
	}
}
