package contraband;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.core.Input.Validate;

public class OUMVNLikelihoodOneTrait extends MVNProcessOneTrait {

	final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Sigma^2, the variance of the process.", Validate.REQUIRED);
	final public Input<RealParameter> rootValueInput = new Input<>("rootValue", "Root trait value, or y_0, also the mean of the stationary distribution at the root if one is assumed.", Validate.REQUIRED);
	final public Input<Boolean> eqDistInput = new Input<>("eqDist", "Whether or not to assume equilibrium (stationary) distribution at the root. The mean of that distribution will rootValue", Validate.REQUIRED);
	final public Input<Boolean> useRootMetaDataInput = new Input<>("useRootMetaData", "Whether or not to use root meta data (specified optimum). If set to 'false', root optimum is set to eldest regime (regimes are numbered from the root toward the tips).", Validate.REQUIRED);
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	final public Input<TreeToVCVMat> optimumManagerInput = new Input<>("optimumManager", "color manager object that paints branches with their own optima.", Validate.REQUIRED);
	final public Input<RealParameter> alphaInput = new Input<>("alpha", "Pull toward optimum or optima.", Validate.REQUIRED);
	//final public Input<ColorManager> optimumManagerInput = new Input<>("optimumManager", "color manager object that paints branches with their own optima.", Validate.REQUIRED);
	
	private boolean dirty;
	private boolean eqDistAtRoot;
	private boolean useRootMetaData;
		
	// basic info
    // private List<Node> allLeafNodes;
	private int nSpp; // takes some time to compute, so part of state!
	//private int nOptima; // used in old W mat parameterization (could remove from state technically, but am calling it a few times, so will stay here!)
	
	/* parameters below */
	// private Double alpha; // could remove from state technically, but am calling it a few times, so will stay here!
	// private Double rootValue; // this is y_0, the value we condition on, or assume has a stationary distr
	                          // the value we condition y_0 on is often theta_0 (root optimum), which is a very liberal assumption!
	
	// mean vector
	private TreeToVCVMat optimumManager;
	// private ColorManager optimumManager;
	// private RealMatrix wMat; // needed only for old parameterization with W matrix
	
	// VCV matrix
	// private RealVector thetaVector; // needed only for old parameterization with W matrix
	private RealVector ouMeanVector;
	private RealMatrix ouTMat, ouVCVMat, ouInvVCVMat;
	private LUDecomposition ouVCVMatLUD;
	
	// data
	private OneValueContTraits oneTraitData;
	private RealVector oneTraitDataVector;
	
	// stored stuff
	private RealMatrix storedVCVMat;
	private RealVector storedOUMeanVector;
	private RealVector storedOneTraitDataVector;
	
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
		// if (useRootMetaData) {
		//	 thetaVector = new ArrayRealVector(nOptima+1);
		//	 wMat = new Array2DRowRealMatrix(nSpp, (nOptima+1));
		// }
		// else {
		//	 thetaVector = new ArrayRealVector(nOptima);
		//	 wMat = new Array2DRowRealMatrix(nSpp, nOptima);
		// }
		 
		oneTraitDataVector = new ArrayRealVector(nSpp);
		ouMeanVector = new ArrayRealVector(nSpp);
		ouTMat = new Array2DRowRealMatrix(nSpp, nSpp);
		
		// store	
		storedOneTraitDataVector = new ArrayRealVector(nSpp);
		storedVCVMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
		storedOUMeanVector = new ArrayRealVector(nSpp);
		
		// this instance vars
		populateInstanceVars(true, true, true, true);
	
		// setting parent class instance vars
		populateParentInstanceVars(true, true, true);
	}
	
	private void populateInstanceVars(boolean updatePhyloTMat, boolean updateVCVMat, boolean updateMean, boolean updateData) {
		if (updateMean) { populateMeanVector(); }
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
		if (updateMean) { setProcessMeanVec(ouMeanVector); }
		if (updateVCVMat) {
			setProcessVCVMat(ouVCVMat);
			setProcessInvVCVMat(ouInvVCVMat);
		}
		if (updateData) {
			setProcessOneTraitDataVec(oneTraitDataVector);
		}
	}
	
	private void populateOUTMatrix(boolean useEqDistBool) {
		double alpha = alphaInput.get().getValue().doubleValue();
		OUUtils.computeOUTMatOneTrait(nSpp, alpha, getPhyloTMatDouble(), ouTMat, useEqDistBool);
	}
	
	@Override
	protected void populateMeanVector() {
		Double rootValue = rootValueInput.get().getValue();		
		// Integer[] thetaAssignments = optimumManager.getColorAssignments();
		// Double[] thetas = optimumManager.getColorValues();
		
		double alpha = alphaInput.get().getValue().doubleValue();
		
		OUUtils.populateOUMeanVector(alpha, rootValue, getRootNode(), getRootNode(), getAllLeafNodes(), optimumManager.getClockModel(), ouMeanVector, useRootMetaData);
		// System.out.println("Mean vector:");
		// GeneralUtils.displayRealVector(ouMeanVector);	

		/*
		 * Old implementation with W matrix
		 */
//			int i = 0;
//			if (useRootMetaData) {
//				thetaVector.setEntry(0, rootValue);
//				i++;
//			}
//			for (Double aTheta: thetas) {
//				thetaVector.setEntry(i, aTheta);
//				i++;
//			} 
//
//			resetRealMatrix(wMat);
//			OUUtils.computeWMatOneTrait(thetaAssignments, getRootNode(), getAllLeafNodes(), nSpp, nOptima, alpha, wMat, useRootMetaData);
//			
//			System.out.println("W matrix:");
//			GeneralUtils.displayRealMatrix(wMat);
//			System.out.println("Theta vector:");
//			GeneralUtils.displayRealVector(thetaVector);
//			System.out.println("Mean vector:");
//			GeneralUtils.displayRealVector(wMat.operate(thetaVector));
//			
//			ouMeanVector = wMat.operate(thetaVector);	
	
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
		oneTraitData = oneTraitInput.get();
		
		int i = 0;
		for (Double traitValue: oneTraitData.getTraitValues(0, getSpNamesInPhyloTMatOrder())) {
			oneTraitDataVector.setEntry(i, traitValue);
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
			storedOUMeanVector.setEntry(i, ouMeanVector.getEntry(i)); // for JIVE
			storedOneTraitDataVector.setEntry(i, oneTraitDataVector.getEntry(i));

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
		
		realVecTmp = ouMeanVector;
		ouMeanVector = storedOUMeanVector;
		storedOUMeanVector = realVecTmp;
		
		realVecTmp = oneTraitDataVector;
		oneTraitDataVector = storedOneTraitDataVector;
		storedOneTraitDataVector = realVecTmp;

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
