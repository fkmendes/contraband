package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import contraband.MVNUtils;

public class BMOneTraitAnalyticalTest {
	
	final static double EPSILON = 1e-4;
	private double resLk; 
	
	@Before
	public void setUP() throws Exception {
		
		// var, mean, n
		int n = 4;
		double var = 1.4822794118;
		double[] meanInput = new double[] { 3.079142, 3.079142, 3.079142, 3.079142 };
		RealVector mean = new ArrayRealVector(meanInput);
				
		// trait values
		double[] dataInput = new double[] { 4.1, 4.5, 5.9, 0.0 };
		RealVector data = new ArrayRealVector(dataInput);
		
		// tree
		double[][] tMatInput = new double[][] {
			{ 3.0, 2.0, 1.0, 0.0 },
			{ 2.0, 3.0, 1.0, 0.0 },
			{ 1.0, 1.0, 3.0, 0.0 },
			{ 0.0, 0.0, 0.0, 3.0 }
			};
			
		RealMatrix tMat = new Array2DRowRealMatrix(tMatInput);
		LUDecomposition tMatLUD = new LUDecomposition(tMat);
		RealMatrix invTMat = tMatLUD.getSolver().getInverse(); // if only variance chances and not tree, we don't have to invert
		RealMatrix invVcvMat = invTMat.scalarMultiply(1/var); // (var*tMat)^-1 = tMat^-1 / var
		double varToNdetTMat = Math.pow(var, n) * tMatLUD.getDeterminant(); // det(var*tMat) = var^n * det(tMat)
		
		resLk = MVNUtils.computeMVNLk(n, var, mean, data, invVcvMat, varToNdetTMat);
		//System.out.println("Likelihood value" + resLk); // 2.4983819506063697E-4
	}
	
	@Test 
	public void againstresLikelihood () {
		Assert.assertEquals(0.0002498382, resLk, EPSILON);
	}


}
