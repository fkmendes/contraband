package test;

import contraband.math.MVNUtils;
import org.junit.Assert;
import org.junit.Test;

public class MVNNormalDensityTest {
	
	final static double EPSILON = 1e-4;
	
	private static double res1, res2, res3;

	// Tests checked with dnorm(x, mean, sd) R function.
	// x: value from -infty to +infty
	// mean: mean of the distribution
	// sd: standard deviation (i.e sqrt(sigma^2))
	// We haven't written an R file for it since it is a basic function
	@Test 
	public void FirstCheck() {
		double[] value = { 0.5, 0.9, -10.5 };
		double[] mean = { 2.3, 0.5, 3 };
		double[] sigsq = { 1, 0.1, 4.2 };

		res1 = MVNUtils.getNormalLk(value[0], mean[0], sigsq[0]);
		res2 = MVNUtils.getNormalLk(value[1], mean[1], sigsq[1]);
		res3 = MVNUtils.getNormalLk(value[2], mean[2], sigsq[2]);
		Assert.assertEquals(0.07895016, res1, EPSILON);	 	// dnorm(0.5, 2.3, 1) 			= 0.07895016
		Assert.assertEquals(0.5668583, res2, EPSILON);		// dnorm(0.9, 0.5, sqrt(0.1)) 	= 0.5668583
		Assert.assertEquals(7.356076e-11, res3, EPSILON);	// dnorm(-10.5, 3, sqrt(4.2)) 	= 7.356076e-11
	}
}