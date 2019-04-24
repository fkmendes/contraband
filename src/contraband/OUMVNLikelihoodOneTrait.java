package contraband;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.sun.javafx.geom.transform.SingularMatrixException;

import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import beast.core.Input.Validate;

public class OUMVNLikelihoodOneTrait extends MVNProcessOneTrait {

	final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Sigma^2, the variance of the process.", Validate.REQUIRED);
	final public Input<RealParameter> rootValueInput = new Input<>("rootValue", "Root trait value, or theta_0, also the mean of the stationary distribution at the root if one is assumed.", Validate.REQUIRED);
	final public Input<Boolean> eqDistInput = new Input<>("eqDist", "Whether or not to assume equilibrium (stationary) distribution at the root. The mean of that distribution will rootValue", Validate.REQUIRED);
	final public Input<Boolean> useRootMetaDataInput = new Input<>("useRootMetaData", "Whether or not to use root meta data (specified optimum). If set to 'false', root optimum is set to eldest regime (regimes are numbered from the root toward the tips).", Validate.REQUIRED);
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	final public Input<ColorManager> optimumManagerInput = new Input<>("optimumManager", "color manager object that paints branches with their own optima.", Validate.REQUIRED);
	final public Input<RealParameter> alphaInput = new Input<>("alpha", "Pull toward optimum or optima.", Validate.REQUIRED);
//	final public Input<RealParameter> thetaInput = new Input<>("theta", "Optimum or optima values, these are the 'colors' of the branches.", Validate.REQUIRED);
//	final public Input<Integer> nOptimaInput = new Input<>("nOptima", "Number of adaptive optima.", Validate.REQUIRED);
	
	private boolean dirty;
	private boolean eqDistAtRoot;
	private boolean useRootMetaData;
	private boolean matrixWasSingularCantInvertBarf;
		
	// basic info
	private List<Node> allLeafNodes;
	private int nSpp; // takes some time to compute, so part of state!
	private int nOptima; // could remove from state technically, but am calling it a few times, so will stay here!
	
	/* parameters below */
	private Double alpha; // could remove from state technically, but am calling it a few times, so will stay here!
	private Double rootValue; // this is y_0, the value we condition on, or assume has a stationary distr
	                          // the value we condition y_0 on is often theta_0 (root optimum), which is a very liberal assumption!
	
	// mean vector
	private ColorManager optimumManager;
	private RealMatrix wMat;
	private Double[] theta;
	private Integer[] thetaAssignments;
	
	// VCV matrix
	private RealVector thetaVector, ouMeanVector;
	private RealMatrix ouTMat, ouVCVMat, ouInvVCVMat;
	private LUDecomposition ouVCVMatLUD;
	
	// data
	private OneValueContTraits oneTraitData;
	private RealVector oneTraitDataVector;
	
	// stored stuff
//	private RealMatrix storedWMat;
	private RealMatrix storedVCVMat;
//	private RealVector storedThetaVector;
//	private Double[] storedTheta;
//	private Integer[] storedThetaAssignments;
	
	private RealVector storedOUMeanVector;
//	private Double storedAlpha;
//	private Double storedRootValue;
	
	@Override
	public void initAndValidate() {	
		
		super.initAndValidate();
		
		// reading in stuff
		allLeafNodes = getRootNode().getAllLeafNodes();
		nSpp = getNSpp();	
		eqDistAtRoot = eqDistInput.get();
		useRootMetaData = useRootMetaDataInput.get();
		optimumManager = optimumManagerInput.get();
		nOptima = optimumManager.getNColors();
		// alpha = alphaInput.get().getValue();
		
		// initializing stuff whose size won't change for now
		if (useRootMetaData) {
			thetaVector = new ArrayRealVector(nOptima+1);
//			storedThetaVector = new ArrayRealVector(nOptima+1);
			wMat = new Array2DRowRealMatrix(nSpp, (nOptima+1));
//			storedWMat = new Array2DRowRealMatrix(nSpp, (nOptima+1));
		}
		else {
			thetaVector = new ArrayRealVector(nOptima);
//			storedThetaVector = new ArrayRealVector(nOptima);
			wMat = new Array2DRowRealMatrix(nSpp, nOptima);
//			storedWMat = new Array2DRowRealMatrix(nSpp, nOptima);
		}
		 
		ouMeanVector = new ArrayRealVector(nSpp);
		ouTMat = new Array2DRowRealMatrix(nSpp, nSpp);
		
		// store
//		storedTheta = new Double[nOptima];
//		storedThetaAssignments = new Integer[optimumManager.getNNodes()];	
		storedVCVMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
		storedOUMeanVector = new ArrayRealVector(nSpp);
		
		// this instance vars
		populateInstanceVars(true, true, true, true);
		populateOneTraitDataVector(); // won't change, so outside populateInstanceVars
		                              // also has to come AFTER setting phyloTMat
		                              // (as this sets the order of species names in T matrix)
		
		// setting parent class instance vars
		populateParentInstanceVars(true, true);
		setProcessOneTraitDataVec(oneTraitDataVector);
	}
	
	private void populateInstanceVars(boolean updatePhyloTMat, boolean updateVCVMat, boolean updateMean, boolean updateAlpha) {
		if (updateAlpha) { alpha = alphaInput.get().getValue(); }
		if (updateMean) { populateMeanVector(); }
		if (updatePhyloTMat) { super.populatePhyloTMatrix(); }
		if (updateVCVMat) {
			populateVCVMatrix();
			populateInvVCVMatrix();
		}
	}
	
	private void populateParentInstanceVars(boolean updateVCVMat, boolean updateMean) {
		// setting parent members
		if (updateMean) { setProcessMeanVec(ouMeanVector); }
		if (updateVCVMat) {
			setProcessVCVMat(ouVCVMat);
			setProcessInvVCVMat(ouInvVCVMat);
		}
	}
	
	private void populateOUTMatrix(boolean useEqDistBool) {
		// alpha has already been dealt in populateInstanceVars
		OUUtils.computeOUTMatOneTrait(nSpp, alpha, getPhyloTMatDouble(), ouTMat, useEqDistBool);
	}
	
	@Override
	protected void populateMeanVector() {
//		if (rootValueInput.get() != null || rootValueInput.isDirty()) {
		rootValue = rootValueInput.get().getValue();
//		}
		
		theta = optimumManager.getColorValues();
		thetaAssignments = optimumManager.getColorAssignments();
		
		int i = 0;
		if (useRootMetaData) {
			thetaVector.setEntry(0, rootValue);
			i++;
		}
		for (Double aTheta: theta) {
			thetaVector.setEntry(i, aTheta);
			i++;
		} 

		resetRealMatrix(wMat);
		OUUtils.computeWMatOneTrait(thetaAssignments, getRootNode(), allLeafNodes, nSpp, nOptima, alpha, wMat, useRootMetaData);

		ouMeanVector = wMat.operate(thetaVector);
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
		oneTraitDataVector = new ArrayRealVector(oneTraitData.getTraitValues(0, getSpNamesInPhyloTMatOrder()));
	}
	
	@Override
	public double calculateLogP() {
		boolean updateAlpha = false;
		boolean updatePhyloTMat = false;
		boolean updateVCVMat = false;
		boolean updateMean = false;
		
		if (alphaInput.isDirty()) { updateAlpha = true; updateVCVMat = true; updateMean = true; }
		if (treeInput.isDirty()) {  updatePhyloTMat = true; updateVCVMat = true; }
		if (sigmasqInput.isDirty() || alphaInput.isDirty()) { updateVCVMat = true; }
		if (rootValueInput.isDirty() || alphaInput.isDirty() || optimumManagerInput.isDirty()) { updateMean = true; }
		
		populateInstanceVars(updatePhyloTMat, updateVCVMat, updateMean, updateAlpha);
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
//		storedRootValue = rootValue;
		
//		System.arraycopy(theta, 0, storedTheta, 0, theta.length);
//		System.arraycopy(thetaAssignments, 0, storedThetaAssignments, 0, thetaAssignments.length);
		
		for (int i=0; i<nSpp; ++i) {
			storedOUMeanVector.setEntry(i, ouMeanVector.getEntry(i));

			for (int j=0; j<nSpp; ++j) {
				storedVCVMat.setEntry(i, j, ouVCVMat.getEntry(i, j));
			}
			
//			for (int j=0; j<wMat.getColumnDimension(); ++j) {
//				storedWMat.setEntry(i, j, wMat.getEntry(i, j));
//			}
		}

//		for (int i=0; i<thetaVector.getDimension(); ++i) {
//			storedThetaVector.setEntry(i, thetaVector.getEntry(i));
//		}
		
		super.store();
	}
	
	@Override
	public void restore() {
		RealVector realVecTmp;
		RealMatrix	realMatTmp;
//		Double[] arrayTmp;
//		Integer[] arrayTmp2;

//		arrayTmp = theta;
//		theta = storedTheta;
//		storedTheta = arrayTmp;
//		
//		arrayTmp2 = thetaAssignments;
//		thetaAssignments = storedThetaAssignments;
//		storedThetaAssignments = arrayTmp2;
//		
//		realMatTmp = wMat;
//		wMat = storedWMat;
//		storedWMat= realMatTmp;
//		
//		realVecTmp = thetaVector;
//		thetaVector = storedThetaVector;
//		storedThetaVector = realVecTmp;
	
		realMatTmp = ouVCVMat;
		ouVCVMat = storedVCVMat;
		storedVCVMat = realMatTmp;
		
		realVecTmp = ouMeanVector;
		ouMeanVector = storedOUMeanVector;
		storedOUMeanVector = realVecTmp;

//		alpha = storedAlpha;
//		rootValue = storedRootValue;
		
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
