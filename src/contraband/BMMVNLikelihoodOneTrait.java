package contraband;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.Input.Validate;

public class BMMVNLikelihoodOneTrait extends MVNProcessOneTrait {

	final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Sigma^2, the variance of the process.", Validate.REQUIRED);
	final public Input<RealParameter> meanInput = new Input<>("mean", "mu, or x_0, the mean of the process (and the values at the root).", Validate.REQUIRED);
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("data", "continuous data values for one trait.", Validate.REQUIRED);
	
	private boolean dirty;
	
	private boolean updatePhyloTMat;
	private boolean updateVCVMat;
	private boolean updateMean;
	
	private double sigmasq;

	private Double bmSingleMeanValue;
	private RealVector bmMeanVector;
	private RealMatrix bmVCVMat, bmInvVCVMat;
	private LUDecomposition bmVCVMatLUD;
	
	private OneValueContTraits oneTraitData;
	private RealVector oneTraitDataVector;
	
	// stored stuff
	private double storedSigmasq;
	private Double storedBMSingleMeanValue;
	private RealVector storedBMMeanVector;
	private RealMatrix storedBMVCVMat, storedBMInvVCVMat;
	
	@Override
	public void initAndValidate() {	
		
		super.initAndValidate();
		
		bmMeanVector = new ArrayRealVector(getNSpp());
		
		sigmasq = sigmasqInput.get().getValue();
		bmSingleMeanValue = meanInput.get().getValue();
		oneTraitData = oneTraitInput.get();
		
		// this instance vars
		populateOneTraitDataVector(); // won't change, so outside populateInstanceVars
		populateInstanceVars(true, true, true);
		
		// setting parent class instance vars
		populateParentInstanceVars(); 
	}
	
	private void populateInstanceVars(boolean updatePhyloTMat, boolean updateVCVMat, boolean updateMean) {
		if (updatePhyloTMat) { super.populatePhyloTMatrix(); }
		if (updateMean) { populateMeanVector(); }
		if (updateVCVMat) {
			populateVCVMatrix();
			populateInvVCVMatrix();
		}
	}
	
	private void populateParentInstanceVars() {
		// setting parent members
		setProcessMeanVec(bmMeanVector);
		setProcessVCVMat(bmVCVMat);
		setProcessInvVCVMat(bmInvVCVMat);
		setProcessOneTraitDataVec(oneTraitDataVector);
	}
	
	@Override
	protected void populateMeanVector() {
		bmMeanVector.set(bmSingleMeanValue);
	}
	
	@Override
	protected void populateVCVMatrix() {
		bmVCVMat = getPhyloTMat().scalarMultiply(sigmasq);
	}
	
	@Override
	protected void populateInvVCVMatrix() {
		bmVCVMatLUD = new LUDecomposition(bmVCVMat);
		bmInvVCVMat = bmVCVMatLUD.getSolver().getInverse();
	}
	
	@Override
	protected void populateOneTraitDataVector() {
		oneTraitDataVector = new ArrayRealVector(oneTraitData.getTraitValues(0));
	}
	
	@Override
	protected void populateLogP() {
		updatePhyloTMat = false;
		updateVCVMat = false;
		updateMean = false;
		
		if (treeInput.isDirty()) {  updatePhyloTMat = true; updateVCVMat = true; }
		if (sigmasqInput.isDirty()) { updateVCVMat = true; }
		if (meanInput.isDirty()) { updateMean = true; }
			
		populateInstanceVars(updatePhyloTMat, updateVCVMat, updateMean);
		
		super.populateLogP();
	}
	
	@Override
	public boolean requiresRecalculation() {
		dirty = super.requiresRecalculation(); // if tree changed, dirty=true from parent, else dirty=false
		
		if (sigmasqInput.isDirty() || meanInput.isDirty()) {
			dirty = true;
		}
		
		return dirty;
	}
	
	@Override
	public void store() {
		for (int i=0; i<getNSpp(); ++i) {
			storedBMMeanVector.setEntry(i, bmMeanVector.getEntry(i));
			
			for (int j=0; j<nSpp; ++j) {
				storedBMVCVMat.setEntry(i, j, bmVCVMat.getEntry(i, j));
				storedBMInvVCVMat.setEntry(i, j, bmInvVCVMat.getEntry(i, j));
			}
		}
		
		// not storing LUDecomposition... hope this is not a problem later.
		super.store();
	}
	
	@Override
	public void restore() {
		RealVector realVecTmp;
		RealMatrix realMatTmp;
		
		realVecTmp = bmMeanVector;
		bmMeanVector = storedBMMeanVector;
		storedBMMeanVector = realVecTmp;
		
		// finish this
		
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
