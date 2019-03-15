import java.util.ArrayList;
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
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import contraband.MVNOneTraitLikelihood;
import contraband.MVNUtils;

public class BMMVNLikelihood extends MVNOneTraitLikelihood {

	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Sigma^2, the variance of the process.", Validate.REQUIRED);
	final public Input<RealParameter> meanInput = new Input<>("mean", "mu, or x_0, the mean of the process (and the values at the root).", Validate.REQUIRED);
	
	private Tree tree;
	private double sigmasq;
	private Double[] meanVector;
	private RealVector mean;
	
	// for phylo T matrix
	private double[] nodeToRootPaths;
	private List<Node> leftLeaves = new ArrayList<>();
	private List<Node> rightLeaves = new ArrayList<>();
	
	private double[][] phyloTMatInput;
	RealMatrix phyloTMat;
	LUDecomposition vcvMatLUD;

	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		sigmasq = sigmasqInput.get().getValue();
		meanVector = meanInput.get().getValues();
		
		populateInstanceVars(); 
	}
	
	public void populateInstanceVars() {
		mean = new ArrayRealVector(meanVector);
		MVNUtils.populateTMatrix(tree, nodeToRootPaths, phyloTMatInput, leftLeaves, rightLeaves); // updates phyloTMatInput
		phyloTMat = new Array2DRowRealMatrix(phyloTMatInput);
		
		// setting parent members
		setProcessVCVMat(phyloTMat.scalarMultiply(sigmasq));
		setProcessMeanVec(mean);
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
