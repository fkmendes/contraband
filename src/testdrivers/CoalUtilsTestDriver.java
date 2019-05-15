package testdrivers;

import contraband.CoalUtils;

public class CoalUtilsTestDriver {

	public static void main(String[] args) {
		System.out.println("Tavare's g10,2, (i=10, j=2, t=3) = " + CoalUtils.getTavareGij(10, 2, 3));
		System.out.println("Heled's pn2k, (n=10, k=2, t=3, ne=1) = " + CoalUtils.getHeledGij(10, 2, 3, 1));
	}
}
