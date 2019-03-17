package contraband;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Distribution;

public abstract class MVNProcessOneTrait extends Distribution {
		
	private int nSpp;
	
	// mean vector
	private RealVector processMeanVec;
	
	// VCV matrix
	private RealMatrix vcvMat, invVCVMat;
	private LUDecomposition vcvMatLUDecomposition;
	private double detVCVMat;
	
	// data
	RealVector oneTraitDataVector;
	
	protected void populatePhyloTMatrix() {};
	
	protected void populateMeanVector() {};
	
	protected void populateVCVMatrix() {};
	
	protected void populateInvVCVMatrix() {};
	
	protected void populateOneTraitDataVector() {};
	
	protected void populateLogP() { 
		logP = MVNUtils.getMVNLogLk(nSpp, processMeanVec, oneTraitDataVector, invVCVMat, detVCVMat);
	};
	
	protected void setProcessNSPP(int aNSpp) {
		nSpp = aNSpp;
	}
	
	protected void setProcessVCVMat(RealMatrix aVCVMat) {
		vcvMat = aVCVMat;
		vcvMatLUDecomposition = new LUDecomposition(vcvMat);
		detVCVMat = vcvMatLUDecomposition.getDeterminant();
	};
	
	protected void setProcessInvVCVMat(RealMatrix aInvVCVMat) {
		invVCVMat = aInvVCVMat;
	}
	
	protected void setProcessMeanVec(RealVector aMeanVector) {
		processMeanVec = aMeanVector;
	};
	
	protected void setProcessOneTraitDataVec(RealVector aOneTraitDataVector) {
		oneTraitDataVector = aOneTraitDataVector;
	}
	
	@Override
	public double calculateLogP() {
		populateLogP();
		
		return logP;
	}
}
