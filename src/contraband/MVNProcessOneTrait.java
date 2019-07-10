package contraband;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public abstract class MVNProcessOneTrait extends Distribution {
	
	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<Double> rootEdgeLengthInput = new Input<>("rootEdgeLength", "root edge length.", 0.0, Validate.OPTIONAL);
	
	// TODO: rootEdgeLength is not a parameter, so in the future after I implement the coalescent correction, this should be figured out
	// from the tree (in StarBeast2)
	
	private boolean dirty;
	private boolean matrixWasSingularCantInvertBarf;
	private boolean successiveThetasIncreasing;
	
	// for phylo T matrix
	private Tree tree;
	private int nSpp; // tree has to be traversed, so part of state!
	
	// things below are part of state because I do not want to call new on them every time I compute likelihood
	private double[] nodeToRootPaths;
	private List<Node> leftLeaves;
	private List<Node> rightLeaves;
	private String[] spNamesInPhyloTMatOrder;
	private double[][] phyloTMatDouble; // at some point I should get rid of this and use PhyloTMat straight away
	private RealMatrix phyloTMat, phyloWTMat;
		
	// mean vector
	private RealVector meanVec;
	
	// VCV matrix
	private RealMatrix vcvMat, invVCVMat;
	private LUDecomposition vcvMatLUDecomposition;
	private double detVCVMat;
	
	// ASR stuff 
	private RealVector wVec; // Weight (design) vector for ASR (with multiple traits, it would be a matrix)
	// never changes, so no need to worry about storing
	private RealMatrix ancNodeVCVMat;
	
	// data
	private RealVector oneTraitDataVector;

	// stored stuff	
	private RealVector storedMeanVec; // need for integration with JIVE (unsure why...?)
	private RealMatrix storedPhyloTMat; // needed for tree operators
	private double[][] storedPhyloTMatDouble; // needed for FBD operators
	private RealMatrix storedInvVCVMat; // (below) needed for morphology parameter operators
	private double storedDetVCVMat;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		nSpp = tree.getLeafNodeCount();
		nodeToRootPaths = new double[tree.getNodeCount()];
		leftLeaves = new ArrayList<Node>();
		rightLeaves = new ArrayList<Node>();
		spNamesInPhyloTMatOrder = new String[nSpp];
		phyloTMatDouble = new double[nSpp][nSpp];
		phyloTMat = MatrixUtils.createRealMatrix(phyloTMatDouble);
		phyloWTMat = MatrixUtils.createRealMatrix(tree.getInternalNodeCount(), nSpp);

		// stored stuff
		storedMeanVec = new ArrayRealVector(nSpp);
		storedPhyloTMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
		storedInvVCVMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
		storedPhyloTMatDouble = new double[nSpp][nSpp];
	}
	
	protected void populatePhyloTMatrix() {		
		MVNUtils.populateTMatrix(tree, nodeToRootPaths, phyloTMatDouble, leftLeaves, rightLeaves, spNamesInPhyloTMatOrder); // updates last 3 args
		
		for (int i=0; i<nSpp; ++i) {
			for (int j=0; j<nSpp; ++j) {
				phyloTMat.setEntry(i, j, phyloTMatDouble[i][j]);
			}
		}

		// setting root edge rows/cols in phyloTMat
		Double rootEdgeLength = rootEdgeLengthInput.get();
		if (rootEdgeLength != 0.0) {
			phyloTMat = phyloTMat.scalarAdd(rootEdgeLength);
		}
	};
	
	protected void populateAncNodePhyloTMatrix() {
		MVNUtils.populateAncNodePhyloTMatrix(tree, phyloWTMat);
	}
	
	protected void populateMeanVector() {};
	
	protected void populateVCVMatrix() {};
	
	protected void populateInvVCVMatrix() {};
	
	protected void populateOneTraitDataVector() {};
	
	protected void populateWMatrix() {}; // design matrix used in ASR
	
	protected void populateAncNodeVCVMatrix() {}; // needed for ASR
	
	protected void populateLogP() {
		if (matrixWasSingularCantInvertBarf || !successiveThetasIncreasing || detVCVMat==0.0) {
			myLogP = Double.NEGATIVE_INFINITY;
		}
		else {
			myLogP = MVNUtils.getMVNLogLk(nSpp, meanVec, oneTraitDataVector, invVCVMat, detVCVMat);
		}
	};
	
	// getters
	protected Tree getTree() {
		return tree;
	}
	
	protected int getNSpp() {
		return nSpp;
	}
	
	protected Node getRootNode() {
		return tree.getRoot();
	}
	
	protected List<Node> getAllLeafNodes() {
		return treeInput.get().getRoot().getAllLeafNodes();
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
		return myLogP;
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
	
	protected void setProcessWMat(RealVector aWVec) {
		wVec = aWVec;
	}
	
	protected void setProcessAncNodeVCVMatrix(RealMatrix aAncNodeVCVMat) {
		ancNodeVCVMat = aAncNodeVCVMat;
	}
	
	protected void setMatrixIsSingular(boolean matrixIsSingular) {
		matrixWasSingularCantInvertBarf = matrixIsSingular;
	}
	
	protected void setThetasAreGo(boolean thetasAreGo) {
		successiveThetasIncreasing = thetasAreGo;
	}
	
	// caching
	@Override
	public boolean requiresRecalculation() {
		dirty = false;
		
		if (treeInput.isDirty() || rootEdgeLengthInput.isDirty()) {
			dirty = true;
		}
		
		return dirty;
	}
	
	@Override
	public void store() {	
		for (int i=0; i<nSpp; ++i) {	
			storedMeanVec.setEntry(i, meanVec.getEntry(i));
			
			// necessary for FBD prior
			System.arraycopy(phyloTMatDouble[i], 0, storedPhyloTMatDouble[i], 0, phyloTMatDouble[i].length);
			
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
		RealVector realVecTmp;
		double[][] double2DArrTmp;
		
		realVecTmp = meanVec;
		meanVec = storedMeanVec;
		storedMeanVec = realVecTmp;
		
		realMatTmp = phyloTMat;
		phyloTMat = storedPhyloTMat;
		storedPhyloTMat = realMatTmp;
		
		realMatTmp = invVCVMat;
		invVCVMat = storedInvVCVMat;
		storedInvVCVMat = realMatTmp;
				
		double2DArrTmp = phyloTMatDouble;
		phyloTMatDouble = storedPhyloTMatDouble;
		storedPhyloTMatDouble = double2DArrTmp;
		
		detVCVMat = storedDetVCVMat;
	}
}
