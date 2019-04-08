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

	final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Sigma^2, the variance of the process.", Validate.REQUIRED);
	final public Input<RealParameter> meanInput = new Input<>("mean", "mu, or x_0, the mean of the process (and the values at the root).", Validate.REQUIRED);
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	
	private boolean dirty;
	
	private boolean updatePhyloTMat;
	private boolean updateVCVMat;
	private boolean updateMean;
	
	private int nSpp;
	private double sigmasq;

	private Double bmSingleMeanValue;
	private RealVector bmMeanVector;
	private RealMatrix bmVCVMat, bmInvVCVMat;
	private LUDecomposition bmVCVMatLUD;
	
	private OneValueContTraits oneTraitData;
	private RealVector oneTraitDataVector;
	
	// stored stuff
	private RealVector storedBMMeanVector;
	
	@Override
	public void initAndValidate() {	
		
		super.initAndValidate();
		
		nSpp = getNSpp();
		bmMeanVector = new ArrayRealVector(nSpp);
		oneTraitDataVector = new ArrayRealVector(nSpp);
		storedBMMeanVector = new ArrayRealVector(nSpp);
		
		// this instance vars
		populateInstanceVars(true, true, true);
		// won't change, so outside populateInstanceVars
		                              // also has to come AFTER setting phyloTMat
		                              // (as this sets the order of species names in T matrix)
		
		// setting parent class instance vars
		populateParentInstanceVars(true, true);
	}
	
	private void populateInstanceVars(boolean updatePhyloTMat, boolean updateVCVMat, boolean updateMean) {
		if (updatePhyloTMat) { super.populatePhyloTMatrix(); }
		if (updateMean) { populateMeanVector(); }
		if (updateVCVMat) {
			populateVCVMatrix();
			populateInvVCVMatrix();
		}
		populateOneTraitDataVector();
	}
	
	private void populateParentInstanceVars(boolean updateVCVMat, boolean updateMean) {
		// setting parent members
		if (updateMean) { setProcessMeanVec(bmMeanVector); }
		if (updateVCVMat) {
			setProcessVCVMat(bmVCVMat);
			setProcessInvVCVMat(bmInvVCVMat);
		}
		setProcessOneTraitDataVec(oneTraitDataVector);
	}
	
	@Override
	protected void populateMeanVector() {
		bmSingleMeanValue = meanInput.get().getValue();
		bmMeanVector.set(bmSingleMeanValue);
	}
	
	@Override
	protected void populateVCVMatrix() {
		sigmasq = sigmasqInput.get().getValue();
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
	
	@Override
	public double calculateLogP() {
		updatePhyloTMat = false;
		updateVCVMat = false;
		updateMean = false;

		if (treeInput.isDirty()) {  updatePhyloTMat = true; updateVCVMat = true; }
		if (sigmasqInput.isDirty()) { updateVCVMat = true; }
		if (meanInput.isDirty()) { updateMean = true; }
		
		populateInstanceVars(updatePhyloTMat, updateVCVMat, updateMean);
		populateParentInstanceVars(updateVCVMat, updateMean);
		
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
			storedBMMeanVector.setEntry(i, bmMeanVector.getEntry(i));
		}

		super.store();
	}
	
	@Override
	public void restore() {
		RealVector realVecTmp;
		
//		realVecTmp = oneTraitDataVector;
//		oneTraitDataVector = storedOneTraitDataVector;
//		storedOneTraitDataVector = realVecTmp;
		
		realVecTmp = bmMeanVector;
		bmMeanVector = storedBMMeanVector;
		storedBMMeanVector = realVecTmp;
		
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
