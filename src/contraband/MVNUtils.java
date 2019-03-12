package contraband;

import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class MVNUtils {
	
	/*
	 * Recursive function for producing vcv matrix out of tree
	 */
	public static void fillNodeLeftRightLeaves(Node node, double[] nodeToRootPaths, double[][] tMat, List<Node> leftLeaves, List<Node> rightLeaves) {
		
		int nodeIdx = node.getNr();
		
		if (node.isLeaf()) {
			/*
			 * For comparing it w/ vcv.phylo in R
			 */
			// System.out.println("I am leaf " + node.getID() + " with idx " + node.getNr());
		 
			tMat[nodeIdx][nodeIdx] = nodeToRootPaths[nodeIdx]; // populating diagonal entries (variances)
			return;
		}
		
		if (node.isRoot()) { nodeToRootPaths[nodeIdx] = 0.0;	}
			
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
		
		fillNodeLeftRightLeaves(left, nodeToRootPaths, tMat, leftLeaves, rightLeaves);
		fillNodeLeftRightLeaves(right, nodeToRootPaths, tMat, leftLeaves, rightLeaves);
		
		return;
	}
	
	public static void populateTMatrix(Tree aTree, double[] nodeToRootPaths, double[][] tMat, List<Node> leftLeaves, List<Node> rightLeaves) {
		fillNodeLeftRightLeaves(aTree.getRoot(), nodeToRootPaths, tMat, leftLeaves, rightLeaves);
	}
	
	/*
	 * Main lk computation function
	 */
	
	// author: Joseph Felsenstein
	// title: Maximum likelihood estimation of evolutionary trees from continuous characters
			// Equation 5
	public static double computeMVNLk(int n, double var, RealVector mean, RealVector data, RealMatrix invVcvMat, double varToNdetTMat) {		
		
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
		System.out.println("Printing coming from MVNUtils lik " + likelihood);

		return likelihood;
	}
	
	public static double normalDensity(double x, double mu, double sigma2) {
		
		return (1.0/Math.sqrt(2.0 * Math.PI * sigma2)) * Math.exp( - Math.pow(x - mu, 2)/(2.0 * sigma2));
	}
}
