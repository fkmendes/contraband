package contraband.utils;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;
import java.util.List;

public class MVNUtils {

	// public static Algebra coltAlgebra = new Algebra(); // if we ever care to switch...

	/*
	 * Recursive function for producing T matrix out of tree
	 */
	public static void fillNodeLeftRightLeaves(Node node, double[] nodeToRootPaths, double[][] tMat, List<Node> leftLeaves, List<Node> rightLeaves, String[] spOrderInTMat) {
		
		int nodeIdx = node.getNr();
		
		if (node.isLeaf()) {
			/*
			 * For comparing it w/ vcv.phylo in R
			 */
//			System.out.println("Leaf " + node.getID() + ", nodeIdx=" + nodeIdx); // uncomment to see index of all leaves 
			tMat[nodeIdx][nodeIdx] = nodeToRootPaths[nodeIdx]; // populating diagonal entries (variances)
			spOrderInTMat[nodeIdx] = node.getID();
			return;
		}
		
		if (node.isRoot()) { nodeToRootPaths[nodeIdx] = 0.0; }
			
		// uncomment to see index of internal nodes and which internal node is which
//		System.out.println("Internal node, nodeIdx=" + nodeIdx); // for debugging
//		List<Node> leafNodes = node.getAllLeafNodes();
//		for (Node aNode: leafNodes) {
//			System.out.println("Daughter leaf: " + aNode.getID());
//		}
		
		Node left = node.getChild(0);
		int leftIdx = left.getNr();
		nodeToRootPaths[leftIdx] = nodeToRootPaths[nodeIdx] + left.getLength();
			
		Node right = node.getChild(1);
		int rightIdx = right.getNr();
		nodeToRootPaths[rightIdx] = nodeToRootPaths[nodeIdx] + right.getLength();
		
		if (left.getAllLeafNodes().size() == 0) {
			leftLeaves.clear();
			leftLeaves.add(left);
		} else {
			leftLeaves = left.getAllLeafNodes();
		}
		
		if (right.getAllLeafNodes().size() == 0) {
			rightLeaves.clear();
			rightLeaves.add(right);
		} else {
			rightLeaves = right.getAllLeafNodes();
		}
		
		for (Node aLeftLeaf: leftLeaves) {
			for (Node aRightLeaf: rightLeaves) {
			tMat[aLeftLeaf.getNr()][aRightLeaf.getNr()] = nodeToRootPaths[nodeIdx];
			tMat[aRightLeaf.getNr()][aLeftLeaf.getNr()] = nodeToRootPaths[nodeIdx];
			}
		}
		
		fillNodeLeftRightLeaves(left, nodeToRootPaths, tMat, leftLeaves, rightLeaves, spOrderInTMat);
		fillNodeLeftRightLeaves(right, nodeToRootPaths, tMat, leftLeaves, rightLeaves, spOrderInTMat);
				
		return;
	}
	
	public static void populateTMatrix(Tree aTree, double[] nodeToRootPaths, double[][] tMat, List<Node> leftLeaves, List<Node> rightLeaves, String[] spOrderInTMat) {
		fillNodeLeftRightLeaves(aTree.getRoot(), nodeToRootPaths, tMat, leftLeaves, rightLeaves, spOrderInTMat);
	}
	
	/*
	 * For ASR under univariate BM with single rate
	 * Function populates (in place) the T matrix from the tree for leaves and ancestral nodes
	 */
	public static void populateAncNodePhyloTMatrix(Tree aTree, RealMatrix aPhyloWTMat) {
		double rootHeight = aTree.getRoot().getHeight();
		int nSpp = aTree.getLeafNodeCount();
		
		for (Node leafNode: aTree.getRoot().getAllLeafNodes()) {
			for (Node ancNode: aTree.getInternalNodes()) {
				if (!ancNode.isRoot()) {
					int leafIdx = leafNode.getNr();
					int ancNodeIdx = ancNode.getNr();

					// if ancestral node is immediate parent of leaf, easy, just put in the subtending branch length
					if (ancNode.getAllLeafNodes().contains(leafNode)) {
						aPhyloWTMat.setEntry(ancNodeIdx-nSpp, leafIdx, rootHeight-ancNode.getHeight());
					}
					
					// if ancestral node is NOT immediate parent of leaf, it could still have some covariance with it!
					// we need to recur to the root and see if any node in the anc node's path to the root is a parent
					// of both the anc node and the leaf -- if so, we need to record that subtending branch length
					else {
						Node parentNode = ancNode.getParent();
						
						while (!parentNode.isRoot()) {
							if (parentNode.getAllLeafNodes().contains(leafNode) && parentNode.getChildren().contains(ancNode)) {
								aPhyloWTMat.setEntry(ancNodeIdx-nSpp, leafIdx, rootHeight-parentNode.getHeight());
								break;
							}
							
							parentNode = parentNode.getParent();
						}
					}
				}
			}
		}
	}
	
	public static void populateNodeRateMatrix(Node node, Double[] nodeRates, RealMatrix rateMatrix) {
		int nodeIdx = node.getNr();
		double nodeRate = nodeRates[nodeIdx];
		
		if (node.isLeaf()) {
			rateMatrix.setEntry(nodeIdx, nodeIdx, nodeRate);
			return;
		}
		
		Node leftChild = node.getChild(0);
		List<Node> leftLeaves = new ArrayList<Node>();
		if (!leftChild.isLeaf()) { leftLeaves = leftChild.getAllLeafNodes(); }
		else { leftLeaves.add(leftChild); }
		populateNodeRateMatrix(leftChild, nodeRates, rateMatrix);
		
		Node rightChild = node.getChild(1);
		List<Node> rightLeaves = new ArrayList<Node>();
		if (!rightChild.isLeaf()) { rightLeaves = rightChild.getAllLeafNodes(); }
		else { rightLeaves.add(rightChild); }
		populateNodeRateMatrix(rightChild, nodeRates, rateMatrix);
		
		if (!node.isRoot()) {		
			for (Node leftLeaf: leftLeaves) {
				int leftLeafIdx = leftLeaf.getNr();
				
				for (Node rightLeaf: rightLeaves) {
					int rightLeafIdx = rightLeaf.getNr();
					rateMatrix.setEntry(leftLeafIdx, rightLeafIdx, nodeRate);
					rateMatrix.setEntry(rightLeafIdx, leftLeafIdx, nodeRate);
				}
			}
		}
		else {
			int nSpp = rateMatrix.getRowDimension();
			for (int i=0; i<(nSpp-1); ++i) {
				rateMatrix.setEntry((nSpp-1), i, nodeRate);
				rateMatrix.setEntry(i, (nSpp-1), nodeRate);
			}
		}
	}
	
	public static void populateRateMatrix(Tree aTree, Double[] nodeRates, RealMatrix rateMatrix) {
		Node rootNode = aTree.getRoot();
		populateNodeRateMatrix(rootNode, nodeRates, rateMatrix);
	}
	
	/*
	 * MVN density function
	 * 
	 * Equation 5 in "Maximum-likelihood estimation of evolutionary trees from continuous characters"
	 * by J. Felsenstein (1973) 
	 * 
	 */
	public static double getMVNLk(int n, RealVector mean, RealVector data, RealMatrix invVcvMat, double varToNdetTMat) {		
		
		/*
		 * This whole thing is the normalizing constant the guarantees
		 * an integral of one (a proper density function)
		 */
		double likelihood = 1.0 / ( Math.pow( (2.0 * Math.PI), n/2.0 ) * Math.pow( varToNdetTMat, 0.5 ) );
		
		/*
		 * Now we multiply by the data stuff
		 */
		likelihood *= Math.exp( 
				invVcvMat
				.preMultiply(data.subtract(mean)).mapMultiply(-0.5)
				.dotProduct(data.subtract(mean))
				);

		return likelihood;
	}
	
	public static double getMVNLogLk (int n, RealVector mean, RealVector data, RealMatrix invVcvMat, double detVcvTMat) {
				
		double loglikelihood = -0.5 * (Math.log(detVcvTMat) + n * Math.log(2.0 * Math.PI));
		
		loglikelihood += invVcvMat.preMultiply(data.subtract(mean)).mapMultiply(-0.5)
				.dotProduct(data.subtract(mean));
		
		// System.out.println("BM log-likelihood=" + loglikelihood);
		
		return loglikelihood;
	}
	
	// public static double getMVNLogLkColt (int n, DoubleMatrix1D mean, DoubleMatrix1D data, DoubleMatrix2D invVcvMat, double detVcvMat) {
	//
	//	 DoubleMatrix1D datMinusMean = data.assign(mean, Functions.minus); // x - mu
	//
	//	 double loglikelihood = -0.5 * (Math.log(detVcvMat) + n * Math.log(2.0 * Math.PI));
	//
	//	 loglikelihood += -0.5 * coltAlgebra.mult(coltAlgebra.mult(invVcvMat, datMinusMean), datMinusMean);
	//
	//	 return loglikelihood;
	// }
	
	/*
	 * One-dimensional, simple normal density
	 */
	public static double getNormalLk(double x, double mu, double sigmaSq) {
		return (1.0/Math.sqrt(2.0 * Math.PI * sigmaSq)) * Math.exp(-Math.pow(x - mu, 2)/(2.0 * sigmaSq));
	}
	
	/*
	 * One-dimensional, simple normal density for n samples (same normal density!), in log space
	 * 
	 * Used in JIVE likelihood
	 */
	public static double getSampleNormalLogLk(Double[] samples, Double mu, Double logSigmaSq) {
		double n = samples.length;
		
		double sumToSubtract = 0.0;
		for (int i=0; i<n; ++i) {
			sumToSubtract += Math.pow(samples[i]-mu, 2);
		}

		return (-((n/2.0) * Math.log(2.0 * Math.PI)) +
			   -((n/2.0) * logSigmaSq) +
			   -(1.0/(2.0 * Math.exp(logSigmaSq))) * sumToSubtract);
	}
	
	/*
	 * One-dimensional, simple normal density for n samples (each sample has its own normal density!), in log space
	 * 
	 * Used in white noise (WN) likelihood
	 */
	public static double getSampleMultipleNormalLogLk(Double[] samples, Double[] mus, Double[] sigmaSqs) {
		double n = samples.length;
		
		double firstTerm = (n/2.0) * Math.log(2.0 * Math.PI);
	
		double secondTerm = 0.0;
		double sumToSubtract = 0.0;
		for (int i=0; i<n; ++i) {
			secondTerm += 0.5 * Math.log(sigmaSqs[i]);
			sumToSubtract += (1.0/(2*sigmaSqs[i])) * Math.pow(samples[i]-mus[i], 2);
		}
		
		double loglikelihood = -firstTerm - secondTerm - sumToSubtract;

		return loglikelihood;
	}
	
	/*
	 * Estimates ancestral state under GLS framework
	 */
	public static double getOneTraitAncState(RealMatrix vcvMat, RealMatrix ancVCVMat, RealVector data, Double[] mu) {
		double estimate = 0.0;
		
		return estimate;
	}
}
