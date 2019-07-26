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
import cern.colt.matrix.linalg.Algebra;

public class MVNShiftProcessOneTrait extends Distribution {
	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	
	// colt
	// protected Algebra alg = new Algebra();
	
	private boolean dirty;
	private boolean matrixWasSingularCantInvertBarf;
	private boolean successiveRatesIncreasing;
	
//	private Tree tree;
	private int nSpp;
	
	// mean vector
	// colt
	// private DoubleMatrix1D meanVec;
	// apache
	private RealVector meanVec;
		
	// VCV matrix
	// colt
	// private DoubleMatrix2D vcvMat, invVCVMat;
	// private double detVCVMat;
	// apache
	private RealMatrix vcvMat, invVCVMat;
	private LUDecomposition vcvMatLUDecomposition;
	private double detVCVMat;
		
	// data
	// colt
	// private DoubleMatrix1D oneTraitDataVec;
	// apache
	private RealVector oneTraitDataVec;

	// stored stuff
	// colt
	// private DoubleMatrix1D storedMeanVec;
	// private DoubleMatrix2D storedInvVCVMat;
	// apache
	private RealVector storedMeanVec; // need for integration with JIVE (unsure why...?)
	private RealMatrix storedInvVCVMat; // (below) needed for morphology parameter operators
	 
	private double storedDetVCVMat;
	
	@Override
	public void initAndValidate() {
		Tree tree = treeInput.get();
		nSpp = tree.getLeafNodeCount();

		// stored stuff
		// colt
		// storedMeanVec = DoubleFactory1D.dense.make(nSpp);
		// storedInvVCVMat = DoubleFactory2D.dense.make(nSpp, nSpp);
		// apache
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
//		else if (!successiveRatesIncreasing) {
//			logP = Double.NEGATIVE_INFINITY;
//		}
		else {
			// colt
			// logP = MVNUtils.getMVNLogLkColt(nSpp, meanVec, oneTraitDataVec, invVCVMat, detVCVMat);
			// apache
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
		// System.out.println("logP=" + logP);
		return logP;
	}
	
	// setters
	// apache
	protected void setProcessVCVMat(RealMatrix aVCVMat) {
	// colt
	// protected void setProcessVCVMat(DoubleMatrix2D aVCVMat) {
		vcvMat = aVCVMat;
		// detVCVMat = alg.det(vcvMat);
		
		// apache
		vcvMatLUDecomposition = new LUDecomposition(vcvMat);
		detVCVMat = vcvMatLUDecomposition.getDeterminant();
		
		// apache commons did not throw a singular matrix exception,
		// but we probably got an aggregation of round-off errors
		// during the computation of this determinant -- which can happen
		// if the matrix is singular
		if (Double.isInfinite(detVCVMat)) {
			matrixWasSingularCantInvertBarf = true; 
		}	
	};
	
	// apache
	protected void setProcessInvVCVMat(RealMatrix aInvVCVMat) {
//	protected void setProcessInvVCVMat(DoubleMatrix2D aInvVCVMat) {
		invVCVMat = aInvVCVMat;
	}
	
	// apache
	protected void setProcessMeanVec(RealVector aMeanVector) {
    // colt
	// protected void setProcessMeanVec(DoubleMatrix1D aMeanVector) {
		meanVec = aMeanVector;
	};
	
	// apache
	protected void setProcessOneTraitDataVec(RealVector aOneTraitDataVector) {
	// colt 
	// protected void setProcessOneTraitDataVec(DoubleMatrix1D aOneTraitDataVector) {
		oneTraitDataVec = aOneTraitDataVector;
	}
	
	protected void setMatrixIsSingular(boolean matrixIsSingular) {
		matrixWasSingularCantInvertBarf = matrixIsSingular;
	}
	
	protected void setRatesAreGo(boolean ratesAreGo) {
		successiveRatesIncreasing = ratesAreGo;
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
			// colt
			// storedMeanVec.set(i, meanVec.get(i));
			// apache
			storedMeanVec.setEntry(i, meanVec.getEntry(i));
			
			for (int j=0; j<nSpp; ++j) {
				// colt
				// storedInvVCVMat.set(i, j, invVCVMat.get(i, j));
				// apache
				storedInvVCVMat.setEntry(i, j, invVCVMat.getEntry(i, j));
			}
		}
		
		storedDetVCVMat = detVCVMat;
	}
	
	@Override
	public void restore() {
		// DoubleMatrix2D realMatTmp;
		// DoubleMatrix1D realVecTmp;
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
