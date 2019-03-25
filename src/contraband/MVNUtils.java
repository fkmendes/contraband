package contraband;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
			// System.out.println("I am leaf " + node.getID() + " with idx " + node.getNr());
		 
			tMat[nodeIdx][nodeIdx] = nodeToRootPaths[nodeIdx]; // populating diagonal entries (variances)
			spOrderInTMat[nodeIdx] = node.getID();
			return;
		}
		
		if (node.isRoot()) { nodeToRootPaths[nodeIdx] = 0.0; }
			
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
	public static double normalDensity(double x, double mu, double sigma2) {
		return (1.0/Math.sqrt(2.0 * Math.PI * sigma2)) * Math.exp(-Math.pow(x - mu, 2)/(2.0 * sigma2));
	}
}
