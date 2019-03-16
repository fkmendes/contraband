package contraband;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public abstract class MVNProcessOneTrait extends Distribution {

	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	
	// for phylo T matrix
	private Tree tree;
	private double[] nodeToRootPaths;
	private List<Node> leftLeaves = new ArrayList<>();
	private List<Node> rightLeaves = new ArrayList<>();
	private double[][] phyloTMatInput;
	RealMatrix phyloTMat;
	
	// for mean vector and VCV matrix
	private int nSpp;
	private RealVector processMeanVec, oneTraitData;
	private RealMatrix vcvMat, invVCVMat;
	private LUDecomposition vcvMatLUDecomposition;
	private double detVCVMat;
	
	// data
	RealVector oneTraitDataVector;
	
	protected void populatePhyloTMatrix() {
		MVNUtils.populateTMatrix(tree, nodeToRootPaths, phyloTMatInput, leftLeaves, rightLeaves); // updates phyloTMatInput
		phyloTMat = new Array2DRowRealMatrix(phyloTMatInput);
	}
	
	protected void populateMeanVector() {};
	
	protected void populateVCVMatrix() {};
	
	protected void populateOneTraitDataVector() {};
	
	protected void populateLogP() { 
		logP = MVNUtils.getMVNLogLk(nSpp, processMeanVec, oneTraitData, invVCVMat, detVCVMat);
	};
	
	// getters
	protected RealMatrix getPhyloTMatrix() {
		return phyloTMat;
	}
	
	// setters
	protected void setProcessVCVMat(RealMatrix aVCVMat) {
		vcvMat = aVCVMat;
		vcvMatLUDecomposition = new LUDecomposition(vcvMat);
		detVCVMat = vcvMatLUDecomposition.getDeterminant();
	};
	
	protected void setProcessMeanVec(RealVector aMeanVector) {
		processMeanVec = aMeanVector;
	};
	
	protected void setProcessOneTraitDataVec(RealVector aOneTraitDataVector) {
		oneTraitDataVector = aOneTraitDataVector;
	}
	
	@Override
	public double calculateLogP() {
		populateLogP();
		
		return logP;
	}
}
