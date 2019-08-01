package contraband;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.hamcrest.core.IsNull;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public abstract class MVNProcessOneTrait extends Distribution {
	
	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<Double> rootEdgeLengthInput = new Input<>("rootEdgeLength", "Root edge length.", 0.0, Validate.OPTIONAL);
	final public Input<Boolean> doCoalCorrectionInput = new Input<>("doCoalCorrection", "Whether or not to perform coalescent correction.", false, Validate.REQUIRED);
	final public Input<CoalCorrection> coalCorrectionInput = new Input<>("coalCorrector", "Calculation node that produces a phylogenetic T matrix from tree and (constant) population sizes of its branches.", Validate.OPTIONAL);
	
	// TODO: rootEdgeLength is not a parameter, so in the future after I implement the coalescent correction, this should be figured out
	// from the tree (in StarBeast2)
	
	private boolean dirty;
	private boolean matrixWasSingularCantInvertBarf;
	
	// for phylo T matrix
	private Tree tree;
	private int nSpp; // tree has to be traversed, so part of state!
	private CoalCorrection coal;
	
	// things below are part of state because I do not want to call new on them every time I compute likelihood
	private double[] nodeToRootPaths;
	private List<Node> leftLeaves;
	private List<Node> rightLeaves;
	private String[] spNamesInPhyloTMatOrder;
	private double[][] phyloTMatDouble; // at some point I should get rid of this and use PhyloTMat straight away
	private RealMatrix phyloTMat, phyloWTMat;
		
	// expectation at tip vector
	private RealVector expAtTipVec;
	
	// VCV matrix
	private RealMatrix vcvMat, invVCVMat;
	private LUDecomposition vcvMatLUDecomposition;
	private double detVCVMat;
	
	// ASR stuff 
	// private RealVector wVec; // Weight (design) vector for ASR (with multiple traits, it would be a matrix)
	// never changes, so no need to worry about storing
	private RealMatrix ancNodeVCVMat;
	
	// data
	private RealVector oneTraitDataVector;

	// stored stuff	
	private RealVector storedOneTraitDataVector = new ArrayRealVector(nSpp);
	private RealVector storedExpAtTipVec; // need for integration with JIVE (unsure why...?)
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
		
		coal = null;
		if (doCoalCorrectionInput.get()) {
			coal = coalCorrectionInput.get();
		}

		// stored stuff
		storedOneTraitDataVector = new ArrayRealVector(nSpp);
		storedExpAtTipVec = new ArrayRealVector(nSpp);
		storedPhyloTMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
		storedInvVCVMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
		storedPhyloTMatDouble = new double[nSpp][nSpp];
	}
	
	protected void populatePhyloTMatrix() {
		Double rootEdgeLength = rootEdgeLengthInput.get();
		
		if (doCoalCorrectionInput.get()) {
			phyloTMatDouble = coal.getCorrectedPhyloTMat(spNamesInPhyloTMatOrder);
		}
		else {
			MVNUtils.populateTMatrix(tree, nodeToRootPaths, phyloTMatDouble, leftLeaves, rightLeaves, spNamesInPhyloTMatOrder); // updates last 3 args
		}
		
		// now populating RealMatrix using double[][]
		for (int i=0; i<nSpp; ++i) {
			for (int j=0; j<nSpp; ++j) {
				phyloTMat.setEntry(i, j, phyloTMatDouble[i][j]);
				
				if (rootEdgeLength != 0.0) {
					phyloTMat.addToEntry(i, j, rootEdgeLength);
				}
			}
		}
	}
	
	protected void populateAncNodePhyloTMatrix() {
		MVNUtils.populateAncNodePhyloTMatrix(tree, phyloWTMat);
	}
	
	protected void populateExpAtTipVector() {};
	
	protected void populateVCVMatrix() {};
	
	protected void populateInvVCVMatrix() {};
	
	protected void populateOneTraitDataVector() {};
	
	protected void populateWMatrix() {}; // design matrix used in ASR
	
	protected void populateAncNodeVCVMatrix() {}; // needed for ASR
	
	protected void populateLogP() {
		// if (matrixWasSingularCantInvertBarf || !successiveThetasIncreasing || detVCVMat==0.0) {
		if (matrixWasSingularCantInvertBarf || detVCVMat==0.0) {
			logP = Double.NEGATIVE_INFINITY;
		}
		else {
			logP = MVNUtils.getMVNLogLk(nSpp, expAtTipVec, oneTraitDataVector, invVCVMat, detVCVMat);
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
		expAtTipVec = aMeanVector;
	};
	
	protected void setProcessOneTraitDataVec(RealVector aOneTraitDataVector) {
		oneTraitDataVector = aOneTraitDataVector;
	}
	
	/*
	 * Used in previous parameterization using W matrix for OU
	 */
//	protected void setProcessWMat(RealVector aWVec) {
//		wVec = aWVec;
//	}
	
	protected void setProcessAncNodeVCVMatrix(RealMatrix aAncNodeVCVMat) {
		ancNodeVCVMat = aAncNodeVCVMat;
	}
	
	protected void setMatrixIsSingular(boolean matrixIsSingular) {
		matrixWasSingularCantInvertBarf = matrixIsSingular;
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
			storedOneTraitDataVector.setEntry(i, oneTraitDataVector.getEntry(i)); // for JIVE
			storedExpAtTipVec.setEntry(i, expAtTipVec.getEntry(i));
			
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
		
		realVecTmp = oneTraitDataVector;
		oneTraitDataVector = storedOneTraitDataVector;
		storedOneTraitDataVector = realVecTmp;

		realVecTmp = expAtTipVec;
		expAtTipVec = storedExpAtTipVec;
		storedExpAtTipVec = realVecTmp;
		
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
