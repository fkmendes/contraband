package contraband.mvnlikelihood;

import java.util.List;
import java.util.Random;

import beast.core.Citation;
import beast.core.Description;
import contraband.utils.OUUtils;
import contraband.clock.TreeToVCVMat;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.Input.Validate;

/**
 * @author Fabio K. Mendes
 */

@Description("Ornstein-Uhlenbeck model (matrix algebra implementation) " +
		"for analyses of continuous trait evolution. Adaptive optima are allowed " +
		"to shift.")
@Citation(value = "Hansen, TF (1997). Stabilizing selection " +
		"and the comparative analysis of adaptation. Evolution 51(5), " +
		"1341-1351.", DOI = " 10.1111/j.1558-5646.1997.tb01457.x",
		year = 1997,
		firstAuthorSurname = "Hansen")
public class OUMVNLikelihoodOneTrait extends MVNProcessOneTrait {

	final public Input<RealParameter> sigmasqInput = new Input<>("sigmaSq", "Sigma^2, the variance of the process.", Validate.REQUIRED);
	final public Input<RealParameter> rootValueInput = new Input<>("rootValue", "Root trait value, or y_0, also the mean of the stationary distribution at the root if one is assumed.", Validate.REQUIRED);
	final public Input<Boolean> eqDistInput = new Input<>("eqDist", "Whether or not to assume equilibrium (stationary) distribution at the root. The mean of that distribution will rootValue", Validate.REQUIRED);
	final public Input<Boolean> useRootMetaDataInput = new Input<>("useRootMetaData", "Whether or not to use root meta data (specified optimum). If set to 'false', root optimum is set to eldest regime (regimes are numbered from the root toward the tips).", Validate.REQUIRED);
	final public Input<TreeToVCVMat> optimumManagerInput = new Input<>("optimumManager", "color manager object that paints branches with their own optima.", Validate.REQUIRED);
	final public Input<RealParameter> alphaInput = new Input<>("alpha", "Pull toward optimum or optima.", Validate.REQUIRED);
	final public Input<RealParameter> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	// final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED); // original implementation (for the above line) with OneValueContTraits data wrapper
	// final public Input<ColorManager> optimumManagerInput = new Input<>("optimumManager", "color manager object that paints branches with their own optima.", Validate.REQUIRED);
	
	private boolean dirty;
	private boolean eqDistAtRoot;
	private boolean useRootMetaData;
		
	// basic info
    // private List<Node> allLeafNodes;
	private int nSpp; // takes some time to compute, so part of state!
	private int nOptima; // used in old W mat parameterization (could remove from state technically, but am calling it a few times, so will stay here!)
	
	/* parameters below */
	// private Double alpha;
	// private Double rootValue; // this is y_0, the value we condition on, or assume has a stationary distr
	                          // the value we condition y_0 on is often theta_0 (root optimum), which is a very liberal assumption!
	
	// mean vector
	private TreeToVCVMat optimumManager;
	// private RealVector thetaVector; // needed only for old parameterization with W matrix
	// private ColorManager optimumManager;
	// private RealMatrix wMat; // needed only for old parameterization with W matrix
	
	// VCV matrix
	// private RealVector thetaVector; // needed only for old parameterization with W matrix
	private RealVector ouMeanVec;
	private RealMatrix ouTMat, ouVCVMat, ouInvVCVMat;
	private LUDecomposition ouVCVMatLUD;
	
	// data
	// private OneValueContTraits oneTraitData; // original implementation with OneValueContTraits data wrapper
	private RealVector oneTraitDataVec;
	
	// stored stuff
	private RealMatrix storedVCVMat;
	private RealVector storedOUMeanVec;
	private RealVector storedOneTraitDataVec;

	@Override
	public void initAndValidate() {	
		
		super.initAndValidate();
		
		// reading in stuff
		nSpp = getNSpp();	
		eqDistAtRoot = eqDistInput.get();
		useRootMetaData = useRootMetaDataInput.get();
		optimumManager = optimumManagerInput.get();
		
		// old parameterization using W mat
		// initializing stuff whose size won't change for now
		// nOptima = optimumManager.getNColors();
//		 if (useRootMetaData) {
//			 thetaVector = new ArrayRealVector(nOptima+1);
//			 wMat = new Array2DRowRealMatrix(nSpp, (nOptima+1));
//		 }
//		 else {
//			 thetaVector = new ArrayRealVector(nOptima);
//			 wMat = new Array2DRowRealMatrix(nSpp, nOptima);
//		 }
		 
		oneTraitDataVec = new ArrayRealVector(nSpp);
		ouMeanVec = new ArrayRealVector(nSpp);
		ouTMat = new Array2DRowRealMatrix(nSpp, nSpp);
		
		// store	
		storedOneTraitDataVec = new ArrayRealVector(nSpp);
		storedVCVMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
		storedOUMeanVec = new ArrayRealVector(nSpp);
		
		// this instance vars
		populateInstanceVars(true, true, true, true);
	
		// setting parent class instance vars
		populateParentInstanceVars(true, true, true);
	}
	
	private void populateInstanceVars(boolean updatePhyloTMat, boolean updateVCVMat, boolean updateMean, boolean updateData) {
		if (updateMean) { populateExpAtTipVector(); }
		if (updatePhyloTMat) { super.populatePhyloTMatrix(); }
		if (updateVCVMat) {
			populateVCVMatrix();
			populateInvVCVMatrix();
		}
		if (updateData) {
			populateOneTraitDataVector();
		}
	}
	
	private void populateParentInstanceVars(boolean updateVCVMat, boolean updateMean,  boolean updateData) {
		// setting parent members
		if (updateMean) { setProcessMeanVec(ouMeanVec); }
		if (updateVCVMat) {
			setProcessLUDec(ouVCVMatLUD);
			setProcessInvVCVMat(ouInvVCVMat);
		}
		if (updateData) {
			setProcessOneTraitDataVec(oneTraitDataVec);
		}
	}
	
	private void populateOUTMatrix(boolean useEqDistBool) {
		double alpha = alphaInput.get().getValue().doubleValue();
		OUUtils.computeOUTMatOneTrait(nSpp, alpha, getPhyloTMatDouble(), ouTMat, useEqDistBool);
	}
	
	@Override
	protected void populateExpAtTipVector() {
		Double rootValue = rootValueInput.get().getValue();

		double alpha = alphaInput.get().getValue().doubleValue();
		
		OUUtils.populateOUMeanVector(alpha, rootValue, getRootNode(), getAllLeafNodes(), optimumManager.getClockModel(), ouMeanVec, useRootMetaData);

		// testing some values by hand
		// RealVector thetaVector = new ArrayRealVector(3);
		// RealMatrix wMat = new Array2DRowRealMatrix(4, 3);
		// Integer[] thetaAssignments = new Integer[] { 1, 1, 1, 0, 1, 1, 0 };

		/*
		 * Alternative implementation with W matrix
		 */
		// Integer[] thetaAssignments = optimumManager.getColorAssignments();
		// Double[] thetas = optimumManager.getColorValues();
		//
		// int i = 0;
		// if (useRootMetaData) {
		//	 thetaVector.setEntry(0, rootValue);
		//	 i++;
		// }
        //
		// Double[] thetas = new Double[] { 0.05779027, 0.19382641 };
        //
		// for (Double aTheta: thetas) {
		//	 thetaVector.setEntry(i, aTheta);
		//	 i++;
		// }
		//
		// resetRealMatrix(wMat);
		// OUUtils.computeWMatOneTrait(thetaAssignments, getRootNode(), getAllLeafNodes(), nSpp, nOptima, alpha, wMat, useRootMetaData);
		//
		// ouMeanVector = wMat.operate(thetaVector);
	}
	
	@Override
	protected void populateVCVMatrix() {
		// alpha has already been dealt in populateInstanceVars
		populateOUTMatrix(eqDistAtRoot);
		Double sigmasq = sigmasqInput.get().getValue();
		ouVCVMat = ouTMat.scalarMultiply(sigmasq);
	}
	
	@Override
	protected void populateInvVCVMatrix() {
		ouVCVMatLUD = new LUDecomposition(ouVCVMat);
		
		try {
			ouInvVCVMat = ouVCVMatLUD.getSolver().getInverse();
			setMatrixIsSingular(false);
		} catch (org.apache.commons.math3.linear.SingularMatrixException e) {
			setMatrixIsSingular(true);
		}
	}
	
	@Override
	protected void populateOneTraitDataVector() {
		/* Original implementation with data wrapper, prior to Parameter having .getValue(aString) */
		// oneTraitData = oneTraitInput.get();
		
		// int i = 0;
		// for (Double traitValue: oneTraitData.getTraitValues(0, getSpNamesInPhyloTMatOrder())) {
		//     oneTraitDataVec.setEntry(i, traitValue);
		//	   ++i;
		// }

		String[] spNamesInPhyloTMatOrder = getSpNamesInPhyloTMatOrder();
		RealParameter oneTraitValues = oneTraitInput.get();
		int i = 0;
		for (String spName: spNamesInPhyloTMatOrder) {
			oneTraitDataVec.setEntry(i, oneTraitValues.getValue(spName));
			++i;
		}
	}
	
	@Override
	public double calculateLogP() {
		boolean updatePhyloTMat = false;
		boolean updateVCVMat = false;
		boolean updateMean = false;
		boolean updateData = false; // for Jive!
		
		if (alphaInput.isDirty()) { updateVCVMat = true; updateMean = true; }
		if (treeInput.isDirty()) {  updatePhyloTMat = true; updateVCVMat = true; }
		if (sigmasqInput.isDirty() || alphaInput.isDirty()) { updateVCVMat = true; }
		if (rootValueInput.isDirty() || alphaInput.isDirty() || optimumManagerInput.isDirty()) { updateMean = true; }
		if (oneTraitInput.isDirty()) { updateData = true; }
		
		populateInstanceVars(updatePhyloTMat, updateVCVMat, updateMean, updateData);
		populateParentInstanceVars(updateVCVMat, updateMean, updateData);
		
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
			storedOUMeanVec.setEntry(i, ouMeanVec.getEntry(i)); // for JIVE
			storedOneTraitDataVec.setEntry(i, oneTraitDataVec.getEntry(i));

			for (int j=0; j<nSpp; ++j) {
				storedVCVMat.setEntry(i, j, ouVCVMat.getEntry(i, j));
			}
		}
		
		super.store();
	}
	
	@Override
	public void restore() {
		RealVector realVecTmp;
		RealMatrix	realMatTmp;
	
		realMatTmp = ouVCVMat;
		ouVCVMat = storedVCVMat;
		storedVCVMat = realMatTmp;
		
		realVecTmp = ouMeanVec;
		ouMeanVec = storedOUMeanVec;
		storedOUMeanVec = realVecTmp;
		
		realVecTmp = oneTraitDataVec;
		oneTraitDataVec = storedOneTraitDataVec;
		storedOneTraitDataVec = realVecTmp;

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
	
	public void resetRealVector(RealVector aRealVector) {
		for (int i=0; i<aRealVector.getDimension(); ++i) {
			aRealVector.setEntry(i, 0.0);
		}
	}
	
	public void resetRealMatrix(RealMatrix aRealMatrix) {
		for (int i=0; i<aRealMatrix.getRowDimension(); ++i) {
			for (int j=0; j<aRealMatrix.getColumnDimension(); ++j) {
				aRealMatrix.setEntry(i, j, 0.0);
			}
		}
	}
}
