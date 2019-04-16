package contraband;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class BMMVNShiftLikelihoodOneTrait extends MVNShiftProcessOneTrait {

	final public Input<RealParameter> meanInput = new Input<>("mean", "mu, or x_0, the mean of the process (and the values at the root).", Validate.REQUIRED);
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	final public Input<ColorManager> rateManagerInput = new Input<>("rateManager", "color manager object that paints branches with their own rates.", Validate.REQUIRED);
		
	private boolean dirty;
	private boolean updateVCVMat;
	private boolean updateMean;
	
	private int nSpp;

	// mean vector
	private Double bmSingleMeanValue;
	private RealVector bmMeanVector;
	
	// VCV matrix
	private ColorManager rateManager;
	private double[][] bmVCVMatDouble;
	private RealMatrix bmVCVMat, bmInvVCVMat;
	private LUDecomposition bmVCVMatLUD;
	
	// data
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
		
		bmVCVMatDouble = new double[nSpp][nSpp];
		bmVCVMat = new Array2DRowRealMatrix(bmVCVMatDouble);
		
		rateManager = rateManagerInput.get();
		
		// this instance vars
		populateInstanceVars(true, true);
				
		// setting parent class instance vars
		populateParentInstanceVars(true, true);
	}
	
	protected void populateInstanceVars(boolean updateVCVMat, boolean updateMean) {
		if (updateMean) { populateMeanVector(); }
		if (updateVCVMat) {
			populateVCVMatrix();
			populateInvVCVMatrix();
		}
		populateOneTraitDataVector();
	}
	
	protected void populateParentInstanceVars(boolean updateVCVMat, boolean updateMean) {
		// setting parent members
		if (updateMean) { setProcessMeanVec(bmMeanVector); }
		if (updateVCVMat) {
			setProcessVCVMat(bmVCVMat);
			setProcessInvVCVMat(bmInvVCVMat);
		}
		setProcessOneTraitDataVec(oneTraitDataVector); // also has to come AFTER setting phyloTMat
        											   // (as this sets the order of species names in T matrix)
	}
	
	@Override
	protected void populateMeanVector() {
		bmSingleMeanValue = meanInput.get().getValue();
		bmMeanVector.set(bmSingleMeanValue);
	}
	
	@Override
	protected void populateVCVMatrix() {
		bmVCVMatDouble = rateManager.getSpColorValuesMatOneTrait();
		
		for (int i=0; i<nSpp; ++i) {
			for (int j=0; j<nSpp; ++j) {
				bmVCVMat.setEntry(i, j, bmVCVMatDouble[i][j]);
			}
		}
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
		for (Double thisTraitValue: oneTraitData.getTraitValues(0, rateManager.getSpNamesInVCVMatOrder())) {
			oneTraitDataVector.setEntry(i, thisTraitValue);
			++i;
		}
	}
	
	@Override
	public double calculateLogP() {
		updateVCVMat = false;
		updateMean = false;

		if (treeInput.isDirty() || rateManagerInput.isDirty()) {  updateVCVMat = true; }
		if (meanInput.isDirty()) { updateMean = true; }
		
		populateInstanceVars(updateVCVMat, updateMean);
		populateParentInstanceVars(updateVCVMat, updateMean);
		
		super.populateLogP();
		
		return getLogP();
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
		
		realVecTmp = bmMeanVector;
		bmMeanVector = storedBMMeanVector;
		storedBMMeanVector = realVecTmp;
		
		super.restore();
	}
	
	@Override
	public boolean requiresRecalculation() {	
		// at the moment, there's no tree caching or anything, so if anything changes
		// gotta recalculate...
		dirty = true;
		return dirty;
	}
}

