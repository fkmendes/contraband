package contraband;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import beast.core.Input.Validate;

public class BMMVNLikelihoodOneTrait extends MVNProcessOneTrait {

	final public Input<TreeParser> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Sigma^2, the variance of the process.", Validate.REQUIRED);
	final public Input<RealParameter> meanInput = new Input<>("mean", "mu, or x_0, the mean of the process (and the values at the root).", Validate.REQUIRED);
	final public Input<OneValueContTraits> oneTraitInput = new Input<>("data", "continuous data values for one trait.", Validate.REQUIRED);
	
	// for phylo T matrix
	int nSpp;
	private TreeParser tree;
	private double[] nodeToRootPaths;
	private List<Node> leftLeaves;
	private List<Node> rightLeaves;
	private double[][] BMPhyloTMatInput;
	private RealMatrix BMPhyloTMat;
	
	private double sigmasq;

	private Double[] BMMeanVectorInput;
	private RealVector BMMeanVector;
	private RealMatrix BMVCVMatrix;
	
	private OneValueContTraits oneTraitData;
	private RealVector oneTraitDataVector;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		nSpp = tree.getLeafNodeCount();
		BMPhyloTMatInput = new double[nSpp][nSpp];
		nodeToRootPaths = new double[tree.getNodeCount()];
		leftLeaves = new ArrayList<>();
		rightLeaves = new ArrayList<>();
		
		sigmasq = sigmasqInput.get().getValue();
		BMMeanVectorInput = meanInput.get().getValues();
		oneTraitData = oneTraitInput.get();
		
		// this instance vars
		populateMeanVector();
		populateVCVMatrix();
		populateOneTraitDataVector();
		
		// setting parent class instance vars
		populateParentInstanceVars(); 
	}
	
	@Override
	protected void populateMeanVector() {
		// later see if I need to move get() from initandvalidate here, and check dirtiness before getting
		// same for other parameters
		
		BMMeanVector = new ArrayRealVector(BMMeanVectorInput);
	}
	
	@Override
	protected void populateVCVMatrix() {
		MVNUtils.populateTMatrix(tree, nodeToRootPaths, BMPhyloTMatInput, leftLeaves, rightLeaves); // updates phyloTMatInput
		BMPhyloTMat = new Array2DRowRealMatrix(BMPhyloTMatInput);
		BMVCVMatrix = BMPhyloTMat.scalarMultiply(sigmasq);
	}
	
	@Override
	protected void populateOneTraitDataVector() {
		oneTraitDataVector = new ArrayRealVector(oneTraitData.getTraitValues(0));
	}
	
	private void populateParentInstanceVars() {
		// setting parent members
		setProcessNSPP(nSpp);
		setProcessMeanVec(BMMeanVector);
		setProcessVCVMat(BMVCVMatrix);
		setProcessOneTraitDataVec(oneTraitDataVector);
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
