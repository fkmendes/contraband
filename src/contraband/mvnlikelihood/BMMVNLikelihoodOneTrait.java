package contraband.mvnlikelihood;
import java.util.List;
import java.util.Random;

import beast.core.Citation;
import beast.core.Description;
import org.apache.commons.math3.linear.*;

import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.Input.Validate;

/**
 * @author Fabio K. Mendes
 */

@Description("Brownian motion model (matrix algebra implementation) " +
		"for analyses of continuous trait evolution.")
@Citation(value = "Felsenstein, J (1985). Maximum-likelihood estimation " +
		"of evolutionary trees from continuous characters. Am. J. Hum. Genet. 25, " +
		"471-492.", DOI = "PMCID: PMC1762641",
		year = 1985,
		firstAuthorSurname = "Felsenstein")
public class BMMVNLikelihoodOneTrait extends MVNProcessOneTrait {

	final public Input<RealParameter> sigmasqInput = new Input<>("sigmaSq", "Sigma^2, the variance of the process.", Validate.REQUIRED); // OPTIONAL because BMMVNShift has ColorManager instead
	final public Input<RealParameter> rootValueInput = new Input<>("rootValue", "rootValue, or y_0, the root value and the expected value at the tips.", Validate.REQUIRED);
	final public Input<RealParameter> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	// final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED); // original implementation (for the above line) with OneValueContTraits data wrapper

	private boolean dirty;
	
	private int nSpp;

	private RealVector bmExpAtTipVector, bmWVector;
	private RealMatrix bmVCVMat, bmInvVCVMat, bmAncNodeVCVMat;
	private LUDecomposition bmVCVMatLUDec;
	
	// private OneValueContTraits oneTraitData; // original implementation with OneValueContTraits data wrapper
	private RealVector oneTraitDataVec;
	
	// stored stuff
	private RealVector storedExpAtTipVec;

	@Override
	public void initAndValidate() {	
		
		super.initAndValidate();
		
		nSpp = getNSpp();
		bmExpAtTipVector = new ArrayRealVector(nSpp);
		bmWVector = new ArrayRealVector(nSpp);
		oneTraitDataVec = new ArrayRealVector(nSpp);

		storedExpAtTipVec = new ArrayRealVector(nSpp);

		// this instance vars
		populateInstanceVars(true, true, true, false, true);
		
		// setting parent class instance vars
		populateParentInstanceVars(true, true, false);
		
		// (ASR) just done once as it does not change for BM
		// populateWMatrix();
		// setProcessWMat(bmWVector);
	}
	
	protected void populateInstanceVars(boolean updatePhyloTMat, boolean updateVCVMat, boolean updateExpAtTip, boolean updateAncNodeVCVMat, boolean updateTipValues) {
		if (updatePhyloTMat) { super.populatePhyloTMatrix(); }
		if (updateExpAtTip) { populateExpAtTipVector(); }
		if (updateVCVMat) {
			populateVCVMatrix();
			populateInvVCVMatrix();
		}
		if (updateTipValues) { populateOneTraitDataVector(); }
		
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
			setProcessLUDec(bmVCVMatLUDec);
			setProcessInvVCVMat(bmInvVCVMat);
		}	
		setProcessOneTraitDataVec(oneTraitDataVec); // also has to come AFTER setting phyloTMat
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
		bmVCVMatLUDec = new LUDecomposition(bmVCVMat);

		try {
			bmInvVCVMat = bmVCVMatLUDec.getSolver().getInverse();
			setMatrixIsSingular(false);
		} catch (org.apache.commons.math3.linear.SingularMatrixException e) {
			setMatrixIsSingular(true);
		}
	}
	
	@Override
	protected void populateOneTraitDataVector() {
		/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
		// oneTraitData = oneTraitInput.get(); // original implementation with OneValueContTraits data wrapper
		// int i = 0;
		// for (Double thisTraitValue: oneTraitData.getTraitValues(0, getSpNamesInPhyloTMatOrder())) {
		// oneTraitDataVector.setEntry(i, thisTraitValue);
		// ++i;
	    // }

		String[] spNamesInPhyloTMatOrder = getSpNamesInPhyloTMatOrder();
		RealParameter oneTraitValues = oneTraitInput.get();
		int i = 0;
		for (String spName: spNamesInPhyloTMatOrder) {
			oneTraitDataVec.setEntry(i, oneTraitValues.getValue(spName));
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
		boolean updateTipValues = false;
		boolean updateAncNodeVCVMat = false; // always false for the moment, only tree logger uses this as true
		
		if (treeInput.isDirty() || (doCoalCorrectionInput.get() && coalCorrectionInput.isDirty())) {  updatePhyloTMat = true; updateVCVMat = true; }
		if (sigmasqInput.isDirty()) { updateVCVMat = true; }
		if (rootValueInput.isDirty()) { updateMean = true; }
		if (oneTraitInput.isDirty()) { updateTipValues = true; }
		
		populateInstanceVars(updatePhyloTMat, updateVCVMat, updateMean, updateAncNodeVCVMat, updateTipValues);
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
			storedExpAtTipVec.setEntry(i, bmExpAtTipVector.getEntry(i));

			// debugging coal correction + MSC
//			for (int j=0; j<nSpp; ++j) {
//				storedBmVCVMat.setEntry(i, j, bmVCVMat.getEntry(i, j));
//				storedBmInvVCVMat.setEntry(i, j, bmInvVCVMat.getEntry(i, j));
//			}
			// storedBmVCVMat.setRowVector(i, bmVCVMat.getRowVector(i));
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
