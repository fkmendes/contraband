package contraband;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class multNormOneTraitAnalytical {
	
	public double computeLk(int n, double var, RealVector mean, RealVector data, RealMatrix invVcvMat, double varToNdetTMat) {		
		
		/*
		 * This whole thing is the normalizing constant the guarantees
		 * an integral of one (a proper density function)
		 */
		double likelihood = 1 / ( Math.pow( (2 * Math.PI), n/2 ) *
				Math.pow( varToNdetTMat, 0.5 ) );
		
		/*
		 * Now we multiply by the data stuff
		 */
		likelihood *= Math.exp( 
				invVcvMat
				.preMultiply(data.subtract(mean)).mapMultiply(-0.5)
				.dotProduct(data.subtract(mean))
				);

		return likelihood;
	}
}
