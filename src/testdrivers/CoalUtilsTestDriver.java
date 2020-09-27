package testdrivers;

import contraband.coalescent.CoalUtils;

/**
 * @author Fabio K. Mendes
 */

public class CoalUtilsTestDriver {

	public static void main(String[] args) {
		System.out.println("Tavare's g10,2, (i=10, j=2, t=3) = " + CoalUtils.getTavareGij(10, 2, 3));
		System.out.println("Heled's pn2k, (n=10, k=2, t=3, pop=1) = " + CoalUtils.getHeledGij(10, 2, 3, 1));
		System.out.println("Mean root height, (k=5, pop=10) = " + CoalUtils.getMeanRootHeight(5, 10));
		System.out.println("Heled's pn2k, (n=2, k=1, t=1.77334814024, d=0.319072) = " + CoalUtils.getHeledGij(2, 1, 1.77334814024, 0.319072));
	}
}
