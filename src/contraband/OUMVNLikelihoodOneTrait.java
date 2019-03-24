package contraband;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.core.Input.Validate;

public class OUMVNLikelihoodOneTrait extends MVNProcessOneTrait {

	final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Sigma^2, the variance of the process.", Validate.REQUIRED);
	final public Input<RealParameter> rootValueInput = new Input<>("rootValue", "Root trait value, or theta_0, also the mean of the stationary distribution at the root if one is assumed.", Validate.OPTIONAL);
	final public Input<Boolean> eqDistInput = new Input<>("eqDist", "Whether or not to assume equilibrium (stationary) distribution at the root. The mean of that distribution will rootValue", Validate.REQUIRED);
	final public Input<RealParameter> alphaInput = new Input<>("alpha", "Pull toward optimum or optima.", Validate.REQUIRED);
	final public Input<RealParameter> thetaInput = new Input<>("theta", "Optimum or optima values, these are the 'colors' of the branches.", Validate.REQUIRED);
	final public Input<Integer> nOptimaInput = new Input<>("nOptima", "Number of adaptive optima.", Validate.REQUIRED);
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	
	private boolean dirty;
	private boolean eqDistAtRoot;
	
	private boolean updatePhyloTMat;
	private boolean updateVCVMat;
	private boolean updateMean;
	
	// basic info
	private Node rootNode;
	private List<Node> allLeafNodes;
	private int nSpp, nOptima;
	
	// parameters
	private double sigmasq;
	private Double rootValue, alpha;
	private Double[] theta;
	
	// to be populated for computation
	private RealVector ouMeanVector;
	private RealMatrix ouTMat, wMat, bmVCVMat, bmInvVCVMat;
	private LUDecomposition bmVCVMatLUD;
	
	// data
	private OneValueContTraits oneTraitData;
	private RealVector oneTraitDataVector;
	
	// stored stuff
	private RealVector storedBMMeanVector;
	
	@Override
	public void initAndValidate() {	
		
		super.initAndValidate();
		
		rootNode = getRootNode();
		allLeafNodes = rootNode.getAllLeafNodes();
		nSpp = getNSpp();
		nOptima = nOptimaInput.get();		
		eqDistAtRoot = eqDistInput.get();
		alpha = alphaInput.get().getValue();
		
		wMat = new Array2DRowRealMatrix(nSpp, nOptima);
		ouMeanVector = new ArrayRealVector(nSpp);
		storedBMMeanVector = new ArrayRealVector(nSpp);
		ouTMat = new Array2DRowRealMatrix(nSpp, nSpp);
		
		// this instance vars
		populateInstanceVars(true, true, true);
		populateOneTraitDataVector(); // won't change, so outside populateInstanceVars
		                              // also has to come AFTER setting phyloTMat
		                              // (as this sets the order of species names in T matrix)
		
		// setting parent class instance vars
		populateParentInstanceVars(true, true);
		setProcessOneTraitDataVec(oneTraitDataVector);
	}
	
	private void populateInstanceVars(boolean updatePhyloTMat, boolean updateVCVMat, boolean updateMean) {
		if (updatePhyloTMat) { super.populatePhyloTMatrix(); }
		if (updateMean) { populateMeanVector(); }
		if (updateVCVMat) {
			populateOUTMatrix(eqDistAtRoot);
			populateVCVMatrix();
			populateInvVCVMatrix();
		}
	}
	
	private void populateParentInstanceVars(boolean updateVCVMat, boolean updateMean) {
		// setting parent members
		if (updateMean) { setProcessMeanVec(ouMeanVector); }
		if (updateVCVMat) {
			setProcessVCVMat(bmVCVMat);
			setProcessInvVCVMat(bmInvVCVMat);
		}
	}
	
	private void populateOUTMatrix(boolean useEqDistBool) {
		OUUtils.computeOUTMatOneTrait(nSpp, alpha, getPhyloTMatDouble(), ouTMat, useEqDistBool);
	}
	
	@Override
	protected void populateMeanVector() {
		if (rootValueInput.get() != null) {
			rootValue = rootValueInput.get().getValue();
		}
		
		theta = thetaInput.get().getValues();
		
		// Do OU stuff here
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
		//TODO: this has to come out with traits in the same order as the species order in the newick format
		//TODO: need to pass in the species order as argument of getTraitValues
		oneTraitDataVector = new ArrayRealVector(oneTraitData.getTraitValues(0, getSpNamesInPhyloTMatOrder()));
		// System.out.println(oneTraitDataVector);
	}
	
	@Override
	public double calculateLogP() {
		updatePhyloTMat = false;
		updateVCVMat = false;
		updateMean = false;
		
		if (treeInput.isDirty()) {  updatePhyloTMat = true; updateVCVMat = true; }
		if (sigmasqInput.isDirty()) { updateVCVMat = true; }
		if (rootValueInput.isDirty() || thetaInput.isDirty()) { updateMean = true; }
		
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
			storedBMMeanVector.setEntry(i, ouMeanVector.getEntry(i));
		}

		super.store();
	}
	
	@Override
	public void restore() {
		RealVector realVecTmp;
		
		realVecTmp = ouMeanVector;
		ouMeanVector = storedBMMeanVector;
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
