package contraband;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class multNormOneTraitAnalytical {

	private double var; // evolutionary rate
	private double[] meanInput;
	private RealVector mean;
	
	private double[] dataInput;
	private RealVector data;
	
	private double[][] vcvMat2DInput; // tree
	private int n; // number of species
	private RealMatrix vcvMat;
	private LUDecomposition vcvMatLU;
	
	// likelihood components
	private RealMatrix varTimesVcvMat;
	private RealVector dataMinusMeanVarTimesVcvMat;
	
	/*
	 * Constructor
	 */
	public multNormOneTraitAnalytical(double var, double[] meanInput, double[] dataInput, double[][] vcvMat2DInput) {
		this.var = var;	
		this.meanInput = meanInput;
		this.mean = new ArrayRealVector(meanInput);
		
		// initializing data
		this.dataInput = dataInput;
		this.data = new ArrayRealVector(dataInput);
		
		// initializing matrix (tree)
		this.vcvMat2DInput = vcvMat2DInput;
		this.n = vcvMat2DInput.length;
		this.vcvMat = new Array2DRowRealMatrix(vcvMat2DInput);
		this.vcvMatLU = new LUDecomposition(vcvMat);
	}
	
	public void computeLk(double var, RealVector mean, RealVector data, RealMatrix vcvMat, int n) {
		double normalizingConstant = 1 / Math.pow( (2 * Math.PI * var), n/2 ) *
				Math.pow( vcvMatLU.getDeterminant(), 0.5 );
			
		varTimesVcvMat = vcvMat.scalarMultiply(var);
		dataMinusMeanVarTimesVcvMat = varTimesVcvMat.preMultiply(data.subtract(mean).mapMultiply(-0.5));
		double dataPart = Math.exp(  );
	}
}
