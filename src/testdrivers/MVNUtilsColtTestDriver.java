package testdrivers;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import contraband.MVNUtils;

public class MVNUtilsColtTestDriver {

	public static void main(String[] args) {
		Algebra alg = new Algebra();
		double[][] matContent = new double[][] { {2.0, 1.0, 0.0}, {1.0, 2.0, 0.0}, {0.0, 0.0, 2.0} };
		
		DoubleMatrix2D vcvMat = DoubleFactory2D.dense.make(matContent).assign(Functions.mult(0.2704762)); // assume sigma^2 = 1
		double detVCVMat = alg.det(vcvMat);
		DoubleMatrix2D invVcvMat = alg.inverse(vcvMat);
		
		double[] datContent = new double[] { 4.1, 4.5, 5.9 };
		DoubleMatrix1D dat = DoubleFactory1D.dense.make(datContent);
		
		double[] meanContent = new double[] { 4.985714, 4.985714, 4.985714 };
		DoubleMatrix1D means = DoubleFactory1D.dense.make(meanContent);
		
		// double loglkColt = MVNUtils.getMVNLogLkColt(3, means, dat, invVcvMat, detVCVMat);
		// System.out.println("MVN log-likelihood using colt = " + loglkColt); // see BMMVNLikelihoodOneTraitTest
		}
}
