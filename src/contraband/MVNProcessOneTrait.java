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
	final public Input<Double> rootEdgeLengthInput = new Input<>("rootEdgeLength", "root edge length.", 0.0, Validate.OPTIONAL);
	
	private boolean dirty;
	
	// for phylo T matrix
	int nSpp;
	private TreeParser tree;
	private double rootEdgeLength;
	private double[] nodeToRootPaths;
	private List<Node> leftLeaves;
	private List<Node> rightLeaves;
	private String[] spNamesInPhyloTMatOrder;
	private double[][] phyloTMatDouble;
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
	private RealMatrix storedPhyloTMat; // needed for tree operators
	private RealMatrix storedInvVCVMat; // (below) needed for morphology parameter operators
	private double storedDetVCVMat;
	
	@Override
	public void initAndValidate() {
		rootEdgeLength = rootEdgeLengthInput.get();
		tree = treeInput.get();
		nSpp = tree.getLeafNodeCount();
		nodeToRootPaths = new double[tree.getNodeCount()];
		leftLeaves = new ArrayList<Node>();
		rightLeaves = new ArrayList<Node>();
		spNamesInPhyloTMatOrder = new String[nSpp];
		phyloTMatDouble = new double[nSpp][nSpp];
		phyloTMat = new Array2DRowRealMatrix(phyloTMatDouble);

		// stored stuff
		storedPhyloTMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
		storedInvVCVMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
	}
	
	protected void populatePhyloTMatrix() {
		MVNUtils.populateTMatrix(tree, nodeToRootPaths, phyloTMatDouble, leftLeaves, rightLeaves, spNamesInPhyloTMatOrder); // updates last 3 args
		
		for (int i=0; i<nSpp; ++i) {
			for (int j=0; j<nSpp; ++j) {
			phyloTMat.setEntry(i, j, phyloTMatDouble[i][j]);
			}
		}

		// setting root edge rows/cols in phyloTMat
		if (rootEdgeLength != 0.0) {
			for(int i=0; i<(nSpp-1); ++i) {
				phyloTMat.setEntry(i, (nSpp-1), rootEdgeLength);
				phyloTMat.setEntry((nSpp-1), i, rootEdgeLength);
			}
		}
	};
	
	protected void populateMeanVector() {};
	
	protected void populateVCVMatrix() {};
	
	protected void populateInvVCVMatrix() {};
	
	protected void populateOneTraitDataVector() {};
	
	protected void populateLogP() { 
		logP = MVNUtils.getMVNLogLk(nSpp, meanVec, oneTraitDataVector, invVCVMat, detVCVMat);
	};
	
	// getters
	protected TreeParser getTree() {
		return tree;
	}
	
	protected int getNSpp() {
		return nSpp;
	}
	
	protected Node getRootNode() {
		return tree.getRoot();
	}
	
	protected String[] getSpNamesInPhyloTMatOrder() {
		return spNamesInPhyloTMatOrder;
	}
	
	protected double[][] getPhyloTMatDouble() {
		return phyloTMatDouble;
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
				storedPhyloTMat.setEntry(i, j, phyloTMat.getEntry(i, j));
				storedInvVCVMat.setEntry(i, j, invVCVMat.getEntry(i, j));
			}
		}
		
		storedDetVCVMat = detVCVMat;
	}
	
	@Override
	public void restore() {
		RealMatrix realMatTmp;
		
		realMatTmp = phyloTMat;
		phyloTMat = storedPhyloTMat;
		storedPhyloTMat = realMatTmp;
		
		realMatTmp = invVCVMat;
		invVCVMat = storedInvVCVMat;
		storedInvVCVMat = realMatTmp;
		
		detVCVMat = storedDetVCVMat;
	}
}
