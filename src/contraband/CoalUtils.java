package contraband;

import java.util.Arrays;

import org.apache.commons.math.util.MathUtils;

public class CoalUtils {
	
	// Expected time to first coalescence as a function of theta (Ne)
	public static double getFirstCoalExpTime(int n, double theta) {
		return theta / MathUtils.binomialCoefficient(n, 2);
	}
	
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
	
	/*
	 * i = nFrom lineages coalesce into j = nTo lineages
	 */
	private static double insideSummation(double tau, int k, int nTo, int nFrom) {
		double inSigma = Math.exp(-k * (k-1) * (tau/2.0)) *
				((2 * k - 1) * Math.pow(-1.0, (k - nTo)) * CoalUtils.geomProgParen(nTo, (k - 1)) * CoalUtils.geomProgBracket(nFrom, k))
						/
						(MathUtils.factorial(nTo)
								* (MathUtils.factorial(k-nTo) * CoalUtils.geomProgParen(nFrom, k)));
		
		return inSigma;
	}
	
	/*
	 * Tavare's gij
	 * 
	 * Prob that i = nFrom lineages coalesce into j = nTo lineages after time tau (in coal. units)
	 */
	public static double getTavareGij(int nFrom, int nTo, double tau) {
		double gij = 0.0;
		
		for (int k=nTo; k<=nFrom; ++k) {
			gij += CoalUtils.insideSummation(tau, k, nTo, nFrom);
		}
		
		return gij;
	}
	
	/*
	 * J. Heled's implementation of Tavare's gij
	 * 
	 * nFrom = n lineages coalesce into nTo = k lineages
	 */
	public static double getHeledGij(int nFrom, int nTo, double t, double ne) {
		
		// nFrom = n
		// nTo = k
		
		// (n-k) rates, backwards
		double[] rates = new double[nFrom - nTo + 1];
		int counter = 0;
		for (int i=nTo; i < (nFrom+1); ++i) {
			rates[counter] = (i*(i-1.0))/2.0;
			counter++;
		}
		
		double[] cis = new double[nFrom - nTo + 1];
		cis[0] = 1.0;
		int cisCounter = 1;
		for (int i=1; i < (nFrom - nTo + 1); ++i) {
			double rate = rates[i];
			
			for (int j=0; j < cisCounter; ++j) {
				cis[j] *= (rate / (rate - rates[j]));
			}

			cis[cisCounter] = -Arrays.stream(cis).sum();
			cisCounter++;
		}
		
		double prob = 0.0;
		for (int i=0; i < rates.length; ++i) {
			prob += (cis[i] * Math.exp((-rates[i]/ne) * t));
		}
		
		return prob;
	}
}
