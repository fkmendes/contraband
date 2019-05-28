package contraband;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class MVNShiftProcessOneTrait extends Distribution {
	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	
	private boolean dirty;
	private boolean matrixWasSingularCantInvertBarf;
	
//	private Tree tree;
	private int nSpp;
	
	// mean vector
	private RealVector meanVec;
		
	// VCV matrix
	private RealMatrix vcvMat, invVCVMat;
	private LUDecomposition vcvMatLUDecomposition;
	private double detVCVMat;
		
	// data
	private RealVector oneTraitDataVec;

	// stored stuff
	private RealVector storedMeanVec; // need for integration with JIVE (unsure why...?)
	private RealMatrix storedInvVCVMat; // (below) needed for morphology parameter operators
	private double storedDetVCVMat;
	
	@Override
	public void initAndValidate() {
		Tree tree = treeInput.get();
		nSpp = tree.getLeafNodeCount();

		// stored stuff
		storedMeanVec = new ArrayRealVector(nSpp);
		storedInvVCVMat = MatrixUtils.createRealMatrix(nSpp, nSpp);
	}

	protected void populateMeanVector() {};
	
	protected void populateVCVMatrix() {};
	
	protected void populateInvVCVMatrix() {};
	
	protected void populateOneTraitDataVector() {};
	
	protected void populateLogP() {
		if (matrixWasSingularCantInvertBarf) {
			logP = Double.NEGATIVE_INFINITY;
		}
		else {
			logP = MVNUtils.getMVNLogLk(nSpp, meanVec, oneTraitDataVec, invVCVMat, detVCVMat);
		}
	};
	
	// getters
	protected Tree getTree() {
		return treeInput.get();
	}
	
	protected int getNSpp() {
		return nSpp;
	}
	
	protected Node getRootNode() {
		Tree tree = treeInput.get();
		return tree.getRoot();
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
		oneTraitDataVec = aOneTraitDataVector;
	}
	
	protected void setMatrixIsSingular(boolean matrixIsSingular) {
		matrixWasSingularCantInvertBarf = matrixIsSingular;
	}
	
	// caching
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
			storedMeanVec.setEntry(i, meanVec.getEntry(i));
			
			for (int j=0; j<nSpp; ++j) {
				storedInvVCVMat.setEntry(i, j, invVCVMat.getEntry(i, j));
			}
		}
		
		storedDetVCVMat = detVCVMat;
	}
	
	@Override
	public void restore() {
		RealMatrix realMatTmp;
		RealVector realVecTmp;
		
		realVecTmp = meanVec;
		meanVec = storedMeanVec;
		storedMeanVec = realVecTmp;
		
		realMatTmp = invVCVMat;
		invVCVMat = storedInvVCVMat;
		storedInvVCVMat = realMatTmp;
		
		detVCVMat = storedDetVCVMat;
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
