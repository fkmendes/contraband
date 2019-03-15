package contraband;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Distribution;

public abstract class MVNOneTraitLikelihood extends Distribution {

	int nSpp;
	double[][] processTMat, processWMat;
	RealVector processMeanVec, oneTraitData;
	RealMatrix vcvMat, invVCVMat;
	LUDecomposition vcvMatLUDecomposition;
	double varToTheNthDet;
	
	public void populatePhyloTMat() {};
	
	public void populateProcessTMat() {};
	
	public void populateProcessInvVCVMat() {};
	
	public void populateProcessMeanVec() {};
	
	// Note: maybe later we want to compute logP straight away inside MVNUtils if underflow happens
	public void populateLogP() { 
		logP = Math.log(MVNUtils.getMVNLk(nSpp, processMeanVec, oneTraitData, invVCVMat, varToTheNthDet));
	};
	
	
}
