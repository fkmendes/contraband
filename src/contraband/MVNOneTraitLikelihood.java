package contraband;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Distribution;

public abstract class MVNOneTraitLikelihood extends Distribution {

	private int nSpp;
	private RealVector processMeanVec, oneTraitData;
	private RealMatrix vcvMat, invVCVMat;
	private LUDecomposition vcvMatLUDecomposition;
	private double detVCVMat;
	
	protected void populateLogP() { 
		logP = MVNUtils.getMVNLogLk(nSpp, processMeanVec, oneTraitData, invVCVMat, detVCVMat);
	};
	
	// setters
	protected void setProcessVCVMat(RealMatrix aVCVMat) {
		vcvMat = aVCVMat;
		vcvMatLUDecomposition = new LUDecomposition(vcvMat);
		detVCVMat = vcvMatLUDecomposition.getDeterminant();
	};
	
	protected void setProcessMeanVec(RealVector aMeanVector) {
		processMeanVec = aMeanVector;
	};
	
	protected void setOneTraitData(RealVector aOneTraitData) {
		oneTraitData = aOneTraitData;
	}
	
	@Override
	public double calculateLogP() {
		populateLogP();
		
		return logP;
	}
}
