package contraband;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class BMMVNShiftLikelihoodOneTrait extends MVNShiftProcessOneTrait {

	final public Input<RealParameter> rootValueInput = new Input<>("rootValue", "rootValue, or y_0, the root value and the expected value at the tips.", Validate.REQUIRED);
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	final public Input<TreeToVCVMat> rateManagerInput = new Input<>("rateManager", "color manager object that paints branches with their own rates.", Validate.REQUIRED);
	
	/*
	 * Deprecated rate class
	 */
	// final public Input<ColorManager> rateManagerInput = new Input<>("rateManager", "color manager object that paints branches with their own rates.", Validate.REQUIRED);
	
	private boolean dirty;
	
	private int nSpp;

	// mean vector
	// colt
	// DoubleMatrix1D bmExpAtTipVector;
	// apache
	private RealVector bmExpAtTipVector;
	
	// VCV matrix
	private TreeToVCVMat rateManager;

	// colt
	// private DoubleMatrix2D bmVCVMat, bmInvVCVMat;
	// apache
	private RealMatrix bmVCVMat, bmInvVCVMat;
	private LUDecomposition bmVCVMatLUD;
	
	// data
	private OneValueContTraits oneTraitData;
	
	// colt
	// private DoubleMatrix1D oneTraitDataVec;
	// apache
	private RealVector oneTraitDataVec;
	
	// stored stuff
	// colt
	// DoubleMatrix1D storedExpAtTipVec;
	// apache
	private RealVector storedExpAtTipVec;
	
	@Override
	public void initAndValidate() {
		
		super.initAndValidate();
		
		nSpp = getNSpp();
		
		// colt 
		// bmExpAtTipVector = DoubleFactory1D.dense.make(nSpp);
		// oneTraitDataVec = DoubleFactory1D.dense.make(nSpp);
		// apache
		bmExpAtTipVector = new ArrayRealVector(nSpp);
		oneTraitDataVec = new ArrayRealVector(nSpp);
				
		// colt
		// bmVCVMat = DoubleFactory2D.dense.make(nSpp, nSpp);
		// apache
		bmVCVMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
		
		rateManager = rateManagerInput.get();
		
		// stored stuff
		// colt
		// storedExpAtTipVec = DoubleFactory1D.dense.make(nSpp);
		// apache
		storedExpAtTipVec = new ArrayRealVector(nSpp);
		
		// this instance vars
		populateInstanceVars(true, true);
				
		// setting parent class instance vars
		populateParentInstanceVars(true, true);
	}
	
	protected void populateInstanceVars(boolean updateVCVMat, boolean updateMean) {
		if (updateMean) { populateExpAtTipVector(); }
		if (updateVCVMat) {
			populateVCVMatrix();
			populateInvVCVMatrix();
		}
		populateOneTraitDataVector();
	}
	
	protected void populateParentInstanceVars(boolean updateVCVMat, boolean updateExpAtTip) {
		// setting parent members
		if (updateExpAtTip) { setProcessMeanVec(bmExpAtTipVector); }
		if (updateVCVMat) {
			setProcessVCVMat(bmVCVMat);
			setProcessInvVCVMat(bmInvVCVMat);
		}
		setProcessOneTraitDataVec(oneTraitDataVec); // also has to come AFTER setting phyloTMat
        											   // (as this sets the order of species names in T matrix)
	}
	
	@Override
	protected void populateExpAtTipVector() {
		Double rootValue = rootValueInput.get().getValue();
		
		// colt
		// for (int i=0; i<nSpp; ++i) {
		//	  bmMeanVector.set(i, bmSingleMeanValue); // repeating the same value nSpp times
		// }
		
		// apache
		bmExpAtTipVector.set(rootValue);
	}
	
	@Override
	protected void populateVCVMatrix() {
		double[][] bmVCVMatDouble = rateManager.getSpColorValuesMatOneTrait();

		for (int i=0; i<nSpp; ++i) {
			for (int j=0; j<nSpp; ++j) {
				// colt
				// bmVCVMat.set(i, j, bmVCVMatDouble[i][j]);
				
				// apache
				bmVCVMat.setEntry(i, j, bmVCVMatDouble[i][j]);
			}
		}
	}
	
	@Override
	protected void populateInvVCVMatrix() {			
		// colt
		// try {
		// 	  bmInvVCVMat = alg.inverse(bmVCVMat);
		//	  setMatrixIsSingular(false);
		// } catch (IllegalArgumentException e) {
		//	  setMatrixIsSingular(true);
		// }
		// apache
		bmVCVMatLUD = new LUDecomposition(bmVCVMat);
					
		try {
			bmInvVCVMat = bmVCVMatLUD.getSolver().getInverse();
			setMatrixIsSingular(false);
		} catch (org.apache.commons.math3.linear.SingularMatrixException e) {
			setMatrixIsSingular(true);
		}
	}
	
	@Override
	protected void populateOneTraitDataVector() {
		oneTraitData = oneTraitInput.get();
		rateManager = rateManagerInput.get();
		
		int i = 0;
		for (Double thisTraitValue: oneTraitData.getTraitValues(0, rateManager.getSpNamesInVCVMatOrder())) {
			// colt 
			// oneTraitDataVec.set(i, thisTraitValue);
			// apache
			oneTraitDataVec.setEntry(i, thisTraitValue);
			++i;
		}
	}
	
	@Override
	public double calculateLogP() {
		boolean updateVCVMat = false;
		boolean updateExpAtTip = false;

		if (treeInput.isDirty() || rateManagerInput.isDirty()) {  updateVCVMat = true; }
		if (rootValueInput.isDirty()) { updateExpAtTip = true; }
		
		populateInstanceVars(updateVCVMat, updateExpAtTip);
		populateParentInstanceVars(updateVCVMat, updateExpAtTip);
		
		super.populateLogP();
		
		return getLogP();
	}
	
	@Override
	public void store() {
		for (int i=0; i<nSpp; ++i) {
			// colt
			// storedExpAtTipVec.set(i, bmExpAtTipVector.get(i));
			// apache
			storedExpAtTipVec.setEntry(i, bmExpAtTipVector.getEntry(i));
		}

		super.store();
	}
	
	@Override
	public void restore() {
		RealVector realVecTmp;

		realVecTmp = bmExpAtTipVector;
		bmExpAtTipVector = storedExpAtTipVec;
		storedExpAtTipVec = realVecTmp;
		
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

