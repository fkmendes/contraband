package contraband;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.Input.Validate;

public class BMMVNLikelihoodOneTrait extends MVNProcessOneTrait {

	final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Sigma^2, the variance of the process.", Validate.REQUIRED); // OPTIONAL because BMMVNShift has ColorManager instead
	final public Input<RealParameter> rootValueInput = new Input<>("rootValue", "rootValue, or y_0, the root value and the expected value at the tips.", Validate.REQUIRED);
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	
	private boolean dirty;
	
	private int nSpp;

	private RealVector bmExpAtTipVector, bmWVector;
	private RealMatrix bmVCVMat, bmInvVCVMat, bmAncNodeVCVMat;
	private LUDecomposition bmVCVMatLUD;
	
	private OneValueContTraits oneTraitData;
	private RealVector oneTraitDataVector;
	
	// stored stuff
	private RealVector storedExpAtTipVector;
	
	@Override
	public void initAndValidate() {	
		
		super.initAndValidate();
		
		nSpp = getNSpp();
		bmExpAtTipVector = new ArrayRealVector(nSpp);
		bmWVector = new ArrayRealVector(nSpp);
		oneTraitDataVector = new ArrayRealVector(nSpp);
		storedExpAtTipVector = new ArrayRealVector(nSpp);
		
		// this instance vars
		populateInstanceVars(true, true, true, false);
		
		// setting parent class instance vars
		populateParentInstanceVars(true, true, false);
		
		// (ASR) just done once as it does not change for BM
		// populateWMatrix();
		// setProcessWMat(bmWVector);
	}
	
	protected void populateInstanceVars(boolean updatePhyloTMat, boolean updateVCVMat, boolean updateExpAtTip, boolean updateAncNodeVCVMat) {
		if (updatePhyloTMat) { super.populatePhyloTMatrix(); }
		if (updateExpAtTip) { populateExpAtTipVector(); }
		if (updateVCVMat) {
			populateVCVMatrix();
			populateInvVCVMatrix();
		}
		populateOneTraitDataVector();
		
		// ASR stuff
		if (updateAncNodeVCVMat) {
			super.populateAncNodePhyloTMatrix();
			populateAncNodeVCVMatrix();
		}
	}
	
	protected void populateParentInstanceVars(boolean updateVCVMat, boolean updateExpAtTip, boolean updateAncNodeVCVMat) {
		// setting parent members
		if (updateExpAtTip) { setProcessMeanVec(bmExpAtTipVector); }
		if (updateVCVMat) {
			setProcessVCVMat(bmVCVMat);
			setProcessInvVCVMat(bmInvVCVMat);
		}	
		setProcessOneTraitDataVec(oneTraitDataVector); // also has to come AFTER setting phyloTMat
        											   // (as this sets the order of species names in T matrix)
		
		// ASR stuff
		if (updateAncNodeVCVMat) {
			setProcessAncNodeVCVMatrix(bmAncNodeVCVMat);
		}
	}
	
	@Override
	protected void populateExpAtTipVector() {
		Double rootValue = rootValueInput.get().getValue();
		bmExpAtTipVector.set(rootValue);
	}
	
	@Override
	protected void populateVCVMatrix() {
		Double sigmasq = sigmasqInput.get().getValue();
		bmVCVMat = getPhyloTMat().scalarMultiply(sigmasq);
	}
	
	@Override
	protected void populateInvVCVMatrix() {
		bmVCVMatLUD = new LUDecomposition(bmVCVMat);
		bmInvVCVMat = bmVCVMatLUD.getSolver().getInverse();
	}
	
	@Override
	protected void populateOneTraitDataVector() {
		oneTraitData = oneTraitInput.get();
		
		int i = 0;
		for (Double thisTraitValue: oneTraitData.getTraitValues(0, getSpNamesInPhyloTMatOrder())) {
			oneTraitDataVector.setEntry(i, thisTraitValue);
			++i;
		}
	}
	
	// ASR stuff
	@Override
	protected void populateWMatrix() {
		bmWVector.set(1.0);
	}
	
	@Override
	protected void populateAncNodeVCVMatrix() {
		// TODO Auto-generated method stub
		;
	}
	
	@Override
	public double calculateLogP() {
		boolean updatePhyloTMat = false;
		boolean updateVCVMat = false;
		boolean updateMean = false;
		boolean updateAncNodeVCVMat = false; // always false, only tree logger uses this as true
		
		if (treeInput.isDirty()) {  updatePhyloTMat = true; updateVCVMat = true; }
		if (sigmasqInput.isDirty()) { updateVCVMat = true; }
		if (rootValueInput.isDirty()) { updateMean = true; }
		
		populateInstanceVars(updatePhyloTMat, updateVCVMat, updateMean, updateAncNodeVCVMat);
		populateParentInstanceVars(updateVCVMat, updateMean, updateAncNodeVCVMat);
		
		super.populateLogP();
		
		return getLogP();
	}
	
	@Override
	public boolean requiresRecalculation() {
//		dirty = super.requiresRecalculation();
//		
//		if (sigmasqInput.isDirty() || meanInput.isDirty()) {
//			dirty = true;
//		}
//
//		return dirty;
		
		// at the moment, there's no tree caching or anything, so if anything changes
		// gotta recalculate...
		dirty = true;
		return dirty;
	}
	
	@Override
	public void store() {
		for (int i=0; i<nSpp; ++i) {
			storedExpAtTipVector.setEntry(i, bmExpAtTipVector.getEntry(i));
		}

		super.store();
	}
	
	@Override
	public void restore() {
		RealVector realVecTmp;
		
		realVecTmp = bmExpAtTipVector;
		bmExpAtTipVector = storedExpAtTipVector;
		storedExpAtTipVector = realVecTmp;
		
		super.restore();
	}
	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
}
