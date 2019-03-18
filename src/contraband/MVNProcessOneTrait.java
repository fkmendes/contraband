package contraband;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.util.TreeParser;

public abstract class MVNProcessOneTrait extends Distribution {
	
	final public Input<TreeParser> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	
	private boolean dirty;
	
	// for phylo T matrix
	int nSpp;
	private TreeParser tree;
	private double[] nodeToRootPaths;
	private List<Node> leftLeaves;
	private List<Node> rightLeaves;
	private double[][] phyloTMatInput;
	private RealMatrix phyloTMat;
		
	// mean vector
	private RealVector meanVec;
	
	// VCV matrix
	private RealMatrix vcvMat, invVCVMat;
	private LUDecomposition vcvMatLUDecomposition;
	private double detVCVMat;
	
	// data
	RealVector oneTraitDataVector;

	// stored stuff
	private RealMatrix storedInvVCVMat;
	private double storedDetVCVMat;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		nSpp = tree.getLeafNodeCount();
		nodeToRootPaths = new double[tree.getNodeCount()];
		leftLeaves = new ArrayList<Node>();
		rightLeaves = new ArrayList<Node>();
		phyloTMatInput = new double[nSpp][nSpp];
		phyloTMat = new Array2DRowRealMatrix(phyloTMatInput);

		// stored stuff
		storedInvVCVMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
	}
	
	protected void populatePhyloTMatrix() {
		MVNUtils.populateTMatrix(tree, nodeToRootPaths, phyloTMatInput, leftLeaves, rightLeaves); // updates last 3 args
		
		for (int i=0; i<nSpp; ++i) {
			for (int j=0; j<nSpp; ++j) {
			phyloTMat.setEntry(i, j, phyloTMatInput[i][j]);
			}
		}

		// phyloTMat = new Array2DRowRealMatrix(phyloTMatInput); // could be used instead of loops
	};
	
	protected void populateMeanVector() {};
	
	protected void populateVCVMatrix() {};
	
	protected void populateInvVCVMatrix() {};
	
	protected void populateOneTraitDataVector() {};
	
	protected void populateLogP() { 
		logP = MVNUtils.getMVNLogLk(nSpp, meanVec, oneTraitDataVector, invVCVMat, detVCVMat);
	};
	
	// getters
	protected int getNSpp() {
		return nSpp;
	}
	
	protected RealMatrix getPhyloTMat() {
		return phyloTMat;
	}
	
	protected double getLogP() {
		return logP;
	}
	
	// setters
	protected void setProcessVCVMat(RealMatrix aVCVMat) {
		vcvMat = aVCVMat;
		vcvMatLUDecomposition = new LUDecomposition(vcvMat);
		detVCVMat = vcvMatLUDecomposition.getDeterminant();
	};
	
	protected void setProcessInvVCVMat(RealMatrix aInvVCVMat) {
		invVCVMat = aInvVCVMat;
	}
	
	protected void setProcessMeanVec(RealVector aMeanVector) {
		meanVec = aMeanVector;
	};
	
	protected void setProcessOneTraitDataVec(RealVector aOneTraitDataVector) {
		oneTraitDataVector = aOneTraitDataVector;
	}
	
	@Override
	public boolean requiresRecalculation() {
		dirty = false;
		
		if (treeInput.isDirty()) {
			dirty = true;
		}
		
		return dirty;
	}
	
	@Override
	public void store() {	
		for (int i=0; i<nSpp; ++i) {		
			for (int j=0; j<nSpp; ++j) {
				storedInvVCVMat.setEntry(i, j, invVCVMat.getEntry(i, j));
			}
		}
		
		storedDetVCVMat = detVCVMat;
	}
	
	@Override
	public void restore() {
		RealMatrix realMatTmp;
		
		realMatTmp = invVCVMat;
		invVCVMat = storedInvVCVMat;
		storedInvVCVMat = realMatTmp;
		
		detVCVMat = storedDetVCVMat;
	}
}
