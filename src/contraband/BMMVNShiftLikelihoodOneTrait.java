package contraband;

//import org.apache.commons.math3.linear.ArrayRealVector;
//import org.apache.commons.math3.linear.LUDecomposition;
//import org.apache.commons.math3.linear.MatrixUtils;
//import org.apache.commons.math3.linear.RealMatrix;
//import org.apache.commons.math3.linear.RealVector;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class BMMVNShiftLikelihoodOneTrait extends MVNShiftProcessOneTrait {

	final public Input<RealParameter> meanInput = new Input<>("mean", "mu, or x_0, the mean of the process (and the values at the root).", Validate.REQUIRED);
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	final public Input<ColorManager> rateManagerInput = new Input<>("rateManager", "color manager object that paints branches with their own rates.", Validate.REQUIRED);
	
	private boolean dirty;
	
	private int nSpp;

	// mean vector
	// colt
	DoubleMatrix1D bmMeanVector;
	// apache
	// private RealVector bmMeanVector;
	
	// VCV matrix
	private ColorManager rateManager;

	// colt
	private DoubleMatrix2D bmVCVMat, bmInvVCVMat;
	// apache
	//	private RealMatrix bmVCVMat, bmInvVCVMat;
	//	private LUDecomposition bmVCVMatLUD;
	
	// data
	private OneValueContTraits oneTraitData;
	
	// colt
	private DoubleMatrix1D oneTraitDataVec;
	// apache
	// private RealVector oneTraitDataVec;
	
	// stored stuff
	// colt
	DoubleMatrix1D storedBMMeanVec;
	// apache
	// private RealVector storedBMMeanVec;
	
	@Override
	public void initAndValidate() {
		
		super.initAndValidate();
		
		nSpp = getNSpp();
		
		// colt 
		bmMeanVector = DoubleFactory1D.dense.make(nSpp);
		oneTraitDataVec = DoubleFactory1D.dense.make(nSpp);
		
		// apache
		// bmMeanVector = new ArrayRealVector(nSpp);
		// oneTraitDataVec = new ArrayRealVector(nSpp);
				
		// colt
		bmVCVMat = DoubleFactory2D.dense.make(nSpp, nSpp);
		// apache
		// bmVCVMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
		
		rateManager = rateManagerInput.get();
		
		// stored stuff
		// colt
		storedBMMeanVec = DoubleFactory1D.dense.make(nSpp);
		// apache
		// storedBMMeanVec = new ArrayRealVector(nSpp);
		
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
		setProcessOneTraitDataVec(oneTraitDataVec); // also has to come AFTER setting phyloTMat
        											   // (as this sets the order of species names in T matrix)
	}
	
	@Override
	protected void populateMeanVector() {
		Double bmSingleMeanValue = meanInput.get().getValue();
		
		// colt
		for (int i=0; i<nSpp; ++i) {
			bmMeanVector.set(i, bmSingleMeanValue); // repeating the same value nSpp times
		}
		
		// apache
		// bmMeanVector.set(bmSingleMeanValue);
	}
	
	@Override
	protected void populateVCVMatrix() {
		double[][] bmVCVMatDouble = rateManager.getSpColorValuesMatOneTrait();

		for (int i=0; i<nSpp; ++i) {
			for (int j=0; j<nSpp; ++j) {
				// colt
				bmVCVMat.set(i, j, bmVCVMatDouble[i][j]);
				
				// apache
				// bmVCVMat.setEntry(i, j, bmVCVMatDouble[i][j]);
			}
		}
	}
	
	@Override
	protected void populateInvVCVMatrix() {
		// colt
		try {
			bmInvVCVMat = alg.inverse(bmVCVMat);
			setMatrixIsSingular(false);
		} catch ( IllegalArgumentException e) {
			setMatrixIsSingular(true);
		}
		
		// apache
//		bmVCVMatLUD = new LUDecomposition(bmVCVMat);
//		
//		try {
//			bmInvVCVMat = bmVCVMatLUD.getSolver().getInverse();
//			setMatrixIsSingular(false);
//		} catch (org.apache.commons.math3.linear.SingularMatrixException e) {
//			setMatrixIsSingular(true);
//		}
	}
	
	@Override
	protected void populateOneTraitDataVector() {
		oneTraitData = oneTraitInput.get();
		
		int i = 0;
		for (Double thisTraitValue: oneTraitData.getTraitValues(0, rateManager.getSpNamesInVCVMatOrder())) {
			// colt 
			oneTraitDataVec.set(i, thisTraitValue);
			// apache
			// oneTraitDataVec.setEntry(i, thisTraitValue);
			++i;
		}
	}
	
	@Override
	public double calculateLogP() {
		boolean updateVCVMat = false;
		boolean updateMean = false;

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
			storedBMMeanVec.set(i, bmMeanVector.get(i));
		}

		super.store();
	}
	
	@Override
	public void restore() {
		DoubleMatrix1D realVecTmp;

		realVecTmp = bmMeanVector;
		bmMeanVector = storedBMMeanVec;
		storedBMMeanVec = realVecTmp;
		
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

