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
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("data", "continuous data values for one trait.", Validate.REQUIRED);
	
	boolean updatePhyloTMat;
	boolean updateVCVMat;
	boolean updateMean;
	
	private double sigmasq;

	private Double BMSingleMeanValue;
	private RealVector BMMeanVector;
	private RealMatrix BMVCVMat, BMInvVCVMat;
	private LUDecomposition BMVCVMatLUD;
	
	private OneValueContTraits oneTraitData;
	private RealVector oneTraitDataVector;
	
	@Override
	public void initAndValidate() {		
		
		super.initAndValidate();
		
		BMMeanVector = new ArrayRealVector(getNSpp());
		
		sigmasq = sigmasqInput.get().getValue();
		BMSingleMeanValue = meanInput.get().getValue();
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
		setProcessMeanVec(BMMeanVector);
		setProcessVCVMat(BMVCVMat);
		setProcessInvVCVMat(BMInvVCVMat);
		setProcessOneTraitDataVec(oneTraitDataVector);
	}
	
	@Override
	protected void populateMeanVector() {
		BMMeanVector.set(BMSingleMeanValue);
	}
	
	@Override
	protected void populateVCVMatrix() {
		BMVCVMat = getPhyloTMat().scalarMultiply(sigmasq);
	}
	
	@Override
	protected void populateInvVCVMatrix() {
		BMVCVMatLUD = new LUDecomposition(BMVCVMat);
		BMInvVCVMat = BMVCVMatLUD.getSolver().getInverse();
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
