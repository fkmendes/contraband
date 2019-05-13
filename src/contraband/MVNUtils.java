package contraband;

import java.util.ArrayList;
import java.util.List;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class MVNUtils {
	
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
	
	public static double getMVNLogLk (int n, RealVector mean, RealVector data, RealMatrix invVcvMat, double detTMat) {
		
		double loglikelihood = -0.5 * (Math.log(detTMat) + n * Math.log(2.0 * Math.PI));
		
		loglikelihood += invVcvMat.preMultiply(data.subtract(mean)).mapMultiply(-0.5)
				.dotProduct(data.subtract(mean));
		
		return loglikelihood;
	}
	
	/*
	 * One-dimensional, simple normal density
	 */
	public static double getNormalLk(double x, double mu, double sigmaSq) {
		return (1.0/Math.sqrt(2.0 * Math.PI * sigmaSq)) * Math.exp(-Math.pow(x - mu, 2)/(2.0 * sigmaSq));
	}
	
	/*
	 * One-dimensional, simple normal density for n samples (same normal density!), in log space
	 */
	public static double getSampleNormalLogLk(double[] samples, Double mu, Double logSigma2) {		
		double n = samples.length;
		
		double sumToSubtract = 0.0;
		for (int i=0; i<n; ++i) {
			sumToSubtract += Math.pow(samples[i]-mu, 2);
		}
		
		return (-((n/2.0) * Math.log(2.0 * Math.PI)) +
			   -((n/2.0) * logSigma2) +
			   -(1.0/(2.0 * Math.exp(logSigma2))) * sumToSubtract);
	}
	
	/*
	 * One-dimensional, simple normal density for n samples (each sample has its own normal density!), in log space
	 * 
	 * Used in white noise (WN) likelihood
	 */
	public static double getSampleMultipleNormalLogLk(Double[] samples, Double[] mus, Double[] logSigmaSqs) {
		double n = samples.length;
		
		double firstTerm = (n/2.0) * Math.log(2.0 * Math.PI);
	
		double secondTerm = 0.0;
		double sumToSubtract = 0.0;
		for (int i=0; i<n; ++i) {
			secondTerm += 0.5 * logSigmaSqs[i];
			sumToSubtract += (1.0/(2*Math.exp(logSigmaSqs[i]))) * Math.pow(samples[i]-mus[i], 2);
		}
		
		return -firstTerm - secondTerm - sumToSubtract;
	}
}
