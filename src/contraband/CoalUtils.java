package contraband;

import java.util.Arrays;

import org.apache.commons.math.util.MathUtils;

public class CoalUtils {
	
	// Expected time to first coalescence as a function of theta (Ne)
	public static double getFirstCoalExpTime(int n, double theta) {
		return theta / MathUtils.binomialCoefficient(n, 2);
	}
	
	/*
	 * Return mean root height when k lineages enter
	 * the root population, which has population size 
	 * pop
	 */
	public static double getMeanRootHeight(int k, double pop) {
		double e = 0.0;
		for (int i=2; i < (k+1); ++i) {
			e += 2.0 * pop / (i * (i-1)); // note the 2.0
		}
		// k * (k-1) is the same as "k take 2"
		
		return e;
	}
	
	/*
	 * Return probability of two lineages coalescing into one
	 * conditioning on a coalescent event within a branch of
	 * length t with population size pop
	 */
	public static double getProbOfTwo2OneCondOnCoal(double t, double pop) {
		return 1 - ( (t/pop) / Math.exp(t/pop)-1 );
	}
	
	/*
	 * Fold distributions of # of starting lineages (k),
	 * one from each daughter population, into the distribution
	 * of the parent (ancestral) population
	 */
	public static double[] foldDistOverNLineages(double[] daughterDist1, double[] daughterDist2) {
		int maxN1 = daughterDist1.length;
		int maxN2 = daughterDist2.length;
		
		double[] foldedParentDist = new double[maxN1 + maxN2];
		int i = 0;
		for (double p1: daughterDist1) {
			
			int j = 0;
			for (double p2: daughterDist2) {
				foldedParentDist[i+j] = p1 * p2;
				j++;
			}
			
			i++;
		}
		
		return foldedParentDist;
	}
	
	/*
	 * Things we need to compute Tavare's Gij's,
	 * where j ancestral lineages evolve into i lineages,
	 * or, conversely, i lineages coalesce into j lineages
	 */
	private static class TavareStuff {
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
					((2 * k - 1) * Math.pow(-1.0, (k - nTo)) * geomProgParen(nTo, (k - 1)) * geomProgBracket(nFrom, k))
							/
							(MathUtils.factorial(nTo)
									* (MathUtils.factorial(k-nTo) * geomProgParen(nFrom, k)));
			
			return inSigma;
		}
	}
	
	/*
	 * Tavare's gij
	 * 
	 * Prob that i = nFrom lineages coalesce into j = nTo lineages after time tau (in coal. units).
	 * Note that, differently from J. Heled's implementation, we have tau in coal. units.!
	 */
	public static double getTavareGij(int nFrom, int nTo, double tau) {		
		double gij = 0.0;
		
		for (int k=nTo; k<=nFrom; ++k) {
			gij += TavareStuff.insideSummation(tau, k, nTo, nFrom);
		}
		
		return gij;
	}
	
	/*
	 * J. Heled's implementation of Tavare's gij
	 * 
	 * Prob that n = nFrom lineages coalesce into k = nTo lineages after time t and given
	 * population size pop.
	 * Note that here we have t in time units and pop as separate parameteres!
	 */
	public static double getHeledGij(int nFrom, int nTo, double t, double pop) {
		
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
			prob += (cis[i] * Math.exp((-rates[i]/pop) * t));
		}
		
		return prob;
	}
}
