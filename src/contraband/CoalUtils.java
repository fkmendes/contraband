package contraband;

import org.apache.commons.math.util.MathUtils;

import com.sun.javafx.css.Combinator;

public class CoalUtils {
	
	// a_(x)
	private static double geomProgParen(int base, int k) {
		if (k == 0) {
			return 1.0;
		}
		
		else {
			int z = 0;
			double prog = 1.0;
			
			while ((base + z) <= (base + k - 1)) {
				prog *= base + z;
				z = z + 1;
			}
			
			return prog;
		}
	}
	
	// a_[x]
	private static double geomProgBracket(int base, int k) {
		if (k == 0) {
			return 1.0;
		}
		
		else {
			int z = 0;
			double prog = 1.0;
			
			while ((base - z) >= (base - k + 1)) {
				prog *= base - z;
				z = z + 1;
			}
			
			return prog;
		}
	}
	
	// inside sigma
	private static double insideSummation(double tau, int k, int j, int i) {
		double inSigma = Math.exp(-k * (k-1) * (tau/2.0)) *
				((2 * k - 1) * Math.pow(-1.0, (k - j)) * CoalUtils.geomProgParen(j, (k - 1)) * CoalUtils.geomProgBracket(i, k))
						/
						(MathUtils.factorial(j)
								* (MathUtils.factorial(k-j) * CoalUtils.geomProgParen(i, k)));
		
		return inSigma;
	}
	
	// Tavare's gij
	public static double getTavareGij(double tau, int j, int i) {
		double gij = 0.0;
		
		for (int k=j; k<=i; ++k) {
			gij += CoalUtils.insideSummation(tau, k, j, i);
		}
		
		return gij;
	}
	
	// Expected time to first coalescence as a function of theta (Ne)
	public static double getFirstCoalExpTime(int n, double theta) {
		return theta / MathUtils.binomialCoefficient(n, 2);
	}
}
