package contraband.utils;

import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import contraband.math.MatrixUtilsContra;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.List;

public class OUUtils {

	public static void computeOUTMatOneTrait(int n, double alpha, double[][] tMat, RealMatrix ouTMat, boolean rootIsRandVar) {	
		double cellValue;
		double divByTwoAlpha = 1.0 / (2.0 * alpha);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				cellValue = divByTwoAlpha * Math.exp(-alpha * (tMat[i][i] + tMat[j][j] - 2.0 * tMat[i][j]));
				
				if (!rootIsRandVar) {
					cellValue *= (1.0 - Math.exp(-2.0 * alpha * tMat[i][j]));
				}
				
				ouTMat.setEntry(i, j, cellValue);  // exponent part
			}
		}	
	}

	public static void computeWMatOneTrait(Integer[] allNodeRegimes, Node rootNode, List<Node> allLeafNodes, int n, int r, double alpha, RealMatrix wMat, boolean useRootMetaData) {
		int rootIndexOffset;
		int rootRegimeIdx;	// Specified in colorAssignment (last position)
			
		if (useRootMetaData) {
			rootIndexOffset = 1; // column dimension of WMat must be r + 1
			rootRegimeIdx = allNodeRegimes[rootNode.getNr()].intValue(); // getting meta data!
		}
		else {
			rootIndexOffset = 0; // column dimension of WMat must be r
			rootRegimeIdx = 0; // not getting metadata, and using the regime that has index 0 (first one)
		}
		 
		for (Node sp: allLeafNodes) {
			int spNr = sp.getNr();	// Will specify the row index

			// Adding the root chunk in the column of the first optimum regime
			double currentSpHeight = sp.getHeight();
			double cellValue = Math.exp(-alpha * (rootNode.getHeight() - currentSpHeight));
			wMat.addToEntry(spNr, rootRegimeIdx, cellValue);
				
			while(!sp.isRoot()) {
				int regimeIdx = allNodeRegimes[sp.getNr()].intValue() + rootIndexOffset; // Will specify the column index
				cellValue = Math.exp(-alpha * (sp.getHeight() - currentSpHeight)) - Math.exp(-alpha * (sp.getParent().getHeight() - currentSpHeight));
				wMat.addToEntry(spNr, regimeIdx, cellValue);

				sp = sp.getParent();
			}
		}
	}
	 
	 /*
	  * Populates OU mean vector by traversing tree only, without pre-computing W matrix
	  * 
	  * Precludes needing to get the number of optima from clock model (since some clock models
	  * do not have discrete rate categories)
	  */
	 public static void populateOUMeanVector(double alpha, double rootValue, Node rootNode, List<Node> allLeafNodes, BranchRateModel clock, RealVector ouMeanVector, boolean useRootMetaData) {
	 // public static void populateOUMeanVector(double alpha, double rootValue, Node aNode, Node rootNode, List<Node> allLeafNodes, BranchRateModel clock, RealVector ouMeanVector, boolean useRootMetaData, double avgSoFar) {
	 	/*
		 * This is a faster, recursive solution, but it only works for ultrametric trees
		 *
		 * It's tricky to generalize this to any tree because the height of a node in BEAST is
		 * always relative to the longest path descending from it, and Hansen's formulation requires
		 * not a branch length(s)=(height1 - height2), but e^{-alpha * height1} - e^{-alpha * height2}.
		 * So note here that for a fossil tip A and a normal tip B below internal node X, X.getHeight()
		 * would have to be larger for B than for A.
		 *
		 * If a tree is non-ultrametric, some of the terms below will be larger than they should be for
		 * extinct tips; e.g., rootNode.getHeight() is a constant, but for extinct tip paths, this
		 * height should actually be smaller (it should be the length going from the root down to that tip
		 * but you can't know that ahead of time during the recursion...)
		 *
		 * Add double avgSoFar to the end parameter list if you want to try to make it work.
		 *
		 */
		// double runningAvg = avgSoFar;

		// if (!aNode.isRoot()) {
		//	 System.out.println("parent height=" + aNode.getParent().getHeight());
		//	 System.out.println("child height=" + aNode.getHeight());
		//	 runningAvg += (Math.exp(-alpha * aNode.getParent().getHeight()) - Math.exp(-alpha * aNode.getHeight())) * clock.getRateForBranch(aNode);
		// }

		// for (Node descNode : aNode.getChildren()) {
		//	 populateOUMeanVector(alpha, rootValue, descNode, rootNode, allLeafNodes, clock, ouMeanVector, useRootMetaData, runningAvg);
		// }

		// if (aNode.isLeaf()) {
		//	 if (useRootMetaData) {
		//		 runningAvg += (Math.exp(-alpha * (rootNode.getHeight() - aNode.getHeight()))) * rootValue;
		//	 }

		//	 ouMeanVector.setEntry(aNode.getNr(), runningAvg);
		// }

		/*
		 * This implementation below follows the notation in Eq. 3 in Hansen 1997.
		 * This is equivalent to Eq. A4 in Butler and King 2004 (there, T = Hansen's tAi).
		 * In Butler & King, if you do the product of the e^{-\alpha T} inside the summation with what's
		 * inside the parenthesis, you end up with:
		 *
		 * (e^{-\alpha * (T - t_i^\gamma)} - e^{-\alpha * (T - t_i^{\gamma -1}})).
		 *
		 * This difference is in fact the same thing as:
		 *
		 *  (e^{-\alpha * tE) - e^{-alpha * tB}).
		 *
		 * Note that in this implementation, I traverse the tree n times, where n is the number of species.
		 * It is slower than the recursive solution above, but works for both ultrametric and non-ultrametric trees
		 *
		 * There's another difference between this implementation and the W * \theta in computeWMatOneTrait;
		 * note that in the latter, \theta_0 = y_0 (there is no independent parameter y_0; y_0 is just the last element in
		 * \theta). We then assign the first optimum regime to the root (i.e., y_0 = \theta_1) when
		 * useRootMetaData = false (i.e., \theta_0 is not used).
		 *
		 * But here, we have a bit more flexibility, and allow the root value y_0 to be its own parameter,
		 * as well as having a separate \theta_0. If useRootMetaData=true, we use y_0 and ignore \theta_0.
		 * Otherwise, \theta_0 can be whatever regime the clock model decides. In the random local clock, \theta_0
		 * propagates down the tree, so it's effectively the "eldest regime" among K regimes. In the RateCategoryClockModel,
		 * \theta_0 will be assigned one of the regimes, not necessarily the first one as in computeWMatOneTrait.
		 * But most importantly, \theta_0 won't be a unique regime in itself (which would be the case of
		 * useRootMetaData=true).
		 *
		 */
		for (Node sp: allLeafNodes) {
			 double cellValue = 0.0;
			 int spNr = sp.getNr();	// Will specify the row index to set in the mean vector

			 double currentSpHeight = sp.getHeight();
			 double tAi = rootNode.getHeight() - currentSpHeight; // tAi is length of path for sp (not simply height of tree)
			 if (useRootMetaData) {
				 cellValue = Math.exp(-alpha * tAi) * rootValue; // use y_0
			 }
			 else {
				 cellValue = Math.exp(-alpha * tAi) * clock.getRateForBranch(rootNode); // use the specified theta for the root
			 }

			// Now going from tip to root, for all tips
			while(!sp.isRoot()) {
				double tE = sp.getHeight() - currentSpHeight;
				double tB = sp.getParent().getHeight() - currentSpHeight; // tB > tE
				cellValue += (Math.exp(-alpha * tE) - Math.exp(-alpha * tB)) * clock.getRateForBranch(sp); // working
				sp = sp.getParent();
			 }

			 ouMeanVector.setEntry(spNr, cellValue);
		 }

		// GeneralUtils.displayRealVector(ouMeanVector);
	 }

	// MultiTrait Covariance matrix functions

	// Pau wrote this, but I should make it faster later...
	public static RealMatrix getExpAlphaMat(int m, RealMatrix eigenDiag, RealMatrix eigenVecMat, double time) {
		RealMatrix res = eigenDiag.createMatrix(m, m);
		
		for (int i = 0; i < m; i++) {
			res.setEntry(i, i, Math.exp( - time * eigenDiag.getEntry(i, i)));
		}
		
		return res.preMultiply(eigenVecMat).multiply(eigenVecMat.transpose());
	}
	
	public static void treeAlphaChunkMat(int m, double[] eigenValues, double cii, double cjj, double cij, boolean fixedRoot, double[][] res) {
		if (fixedRoot) {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < m; j++) {
					res[i][j] = (1 - Math.exp(-(eigenValues[i] + eigenValues[j]) * cij)) / (eigenValues[i] + eigenValues[j]);
				}
			}
		}

		else {
			for (int i = 0; i < m; i++) {
			    for (int j = 0; j < m; j++) {
			    	res[i][j] = (Math.exp( - eigenValues[i] * (cii - cij))
			    			* Math.exp(-eigenValues[j] * (cjj -cij))) / (eigenValues[i] + eigenValues[j]);
			    } 
			}		
		}

	}
	
	public static RealMatrix eigencholChunk(RealMatrix cholSigma2, RealMatrix eigenAlpha) {
		
		return ( eigenAlpha.transpose()
				.multiply(cholSigma2)
				.multiply(cholSigma2.transpose())
				.multiply(eigenAlpha) );
	}
	
	public static RealMatrix multiTraitCovMat2Species(int m, double[][] eigenCholChunk, EigenDecomposition eigenObjectAlpha , double cii, double cjj, double cij, boolean fixedRoot) {
		double[][] res = new double[m][m];
		RealMatrix resMat;
		
		GeneralUtils.displayRealMatrix(eigenObjectAlpha.getD());
		
		if (fixedRoot) {
			
			OUUtils.treeAlphaChunkMat(m, eigenObjectAlpha.getRealEigenvalues(), cii, cjj, cij, fixedRoot, res);
			MatrixUtilsContra.elementWiseProduct(res, eigenCholChunk, res);
			resMat = new Array2DRowRealMatrix(res);

			resMat = OUUtils.getExpAlphaMat(m, eigenObjectAlpha.getD(), eigenObjectAlpha.getV(), cii - cij)
					.multiply(eigenObjectAlpha.getV())
					.multiply(resMat)
					.multiply(eigenObjectAlpha.getVT())
					.multiply(OUUtils.getExpAlphaMat(m, eigenObjectAlpha.getD(), eigenObjectAlpha.getV(), cjj - cij));
		} else {
			
			OUUtils.treeAlphaChunkMat(m, eigenObjectAlpha.getRealEigenvalues(), cii, cjj, cij, false, res);
			MatrixUtilsContra.elementWiseProduct(res, eigenCholChunk, res);
			resMat = new Array2DRowRealMatrix(res);
			resMat = eigenObjectAlpha.getV()
					.multiply(resMat)
					.multiply(eigenObjectAlpha.getVT());
		}
		
		return resMat;
}
	
	
	
	
	// Multitrait Weight Matrix functions
	
	public static RealMatrix multiTraitWeightMatSpVsRegime(int m, Node sp, int k, Node nodeRoot, Double[] regimes, EigenDecomposition EigenAlphaObject, boolean mergeRoot) {
		
		RealMatrix wMat = EigenAlphaObject.getD().createMatrix(m, m); // output matrix
//		Node nodeRoot = Tree.getRoot();
//		List<Node> allLeaf = nodeRoot.getAllLeafNodes();
//		
//		Double[] regimes = new Double[nodeRoot.getNr() + 1]; //It stores the regimes of nodes in the order taken by beast
//		Tree.getMetaData(nodeRoot, regimes, "Regime");	// Populating regimes vector
//		Node sp = allLeaf.get(i); //Selecting the current leaf specified by i

			double spHeight = sp.getHeight();
			
			if(mergeRoot && k == 0) {
				
				RealMatrix cellValue = OUUtils.getExpAlphaMat(m, EigenAlphaObject.getD(), EigenAlphaObject.getV(), nodeRoot.getHeight() - spHeight); 
				wMat = wMat.add(cellValue);
			}

			
			while(!sp.isRoot()) {
				
				int intRegime = regimes[sp.getNr()].intValue();
				System.out.println("regime of species " + intRegime);
				
				if(intRegime == k) {
					
					RealMatrix cellValue = OUUtils.getExpAlphaMat(m, EigenAlphaObject.getD(), EigenAlphaObject.getV(), sp.getHeight() - spHeight)
							.subtract( OUUtils.getExpAlphaMat(m, EigenAlphaObject.getD(), EigenAlphaObject.getV(), sp.getParent().getHeight() - spHeight) );   // Math.exp(-alpha * (sp.getHeight() - currentSpHeight)) - Math.exp(-alpha * (sp.getParent().getHeight() - currentSpHeight));
					wMat = wMat.add(cellValue);
//					GeneralUtils.displayRealMatrix(cellValue);
				}
				sp = sp.getParent();
			}
		
		return wMat;
	}
	
	public static void multiTraitFullWeightMat(int m, Node nodeRoot, Double[] regimes, int r, EigenDecomposition EigenAlphaObject, boolean mergeRoot, double[][] Result) {
		
		List<Node> allLeaf = nodeRoot.getAllLeafNodes();
//		int intRootRegime = regimes[nodeRoot.getNr()].intValue(); //We will set the regime of the root to be always 0
		
		RealMatrix blockWeightMat = EigenAlphaObject.getD().createMatrix(m, m);
		
		if(mergeRoot) {
			
			for(Node sp:allLeaf) {
				
				for(int k = 0; k < r; k++) {
		
					blockWeightMat = multiTraitWeightMatSpVsRegime(m, sp, k, nodeRoot, regimes, EigenAlphaObject, true);
					MatrixUtilsContra.setBlocksInMatrix(blockWeightMat.getData(), m, sp.getNr(), k, Result);
				}
			}
			
		} else {
			
			for(Node sp:allLeaf) {
				
				double spHeight = sp.getHeight();
				RealMatrix cellValue = OUUtils.getExpAlphaMat(m, EigenAlphaObject.getD(), EigenAlphaObject.getV(), nodeRoot.getHeight() - spHeight); 
				MatrixUtilsContra.setBlocksInMatrix(cellValue.getData(), m, sp.getNr(), 0, Result);
				
				for(int k = 0; k < r; k++) {
		
					blockWeightMat = multiTraitWeightMatSpVsRegime(m, sp, k, nodeRoot, regimes, EigenAlphaObject, false);
					MatrixUtilsContra.setBlocksInMatrix(blockWeightMat.getData(), m, sp.getNr(), k + 1, Result);
				}
			}
		}
	}
	
	
//	public static double[][] multiTraitWeightMatSpvsRegime(int m, int k, TreeParser Tree, double[][] phyloMat, EigenDecomposition EigenAlphaObject, boolean mergeRoot) {
//		
//		Node nodeRoot = Tree.getRoot();
//		List<Node> allLeaf = nodeRoot.getAllLeafNodes();
//		
//		Double[] regimes = new Double[nodeRoot.getNr() + 1]; // Will keep the vector with the regimes of the branches subtending each node
//		int rootIndexOffset;
//		int intRootRegime;	// Eldest regime index
//
//		Tree.getMetaData(nodeRoot, regimes, "Regime");
//		
//		double[][] wMat = new double[m][m];
//		if (mergeRoot) { 		// column dimension of WMat must be r
//			
//			rootIndexOffset = 0;
//			intRootRegime = regimes[nodeRoot.getNr()].intValue();
//			
//		} else { 		// column dimension of WMat must be r + 1
//			
//			rootIndexOffset = 1;
//			intRootRegime = 0;	
//		}
//		
//		for (Node sp: allLeaf) {
//				
//			int spNr = sp.getNr();	// Will specify the row index 
//				
//			// Adding the root chunk in the column of the eldest regime
//			double currentSpHeight = sp.getHeight();
//			double[][] cellValue = OUUtils.getExpAlphaMat(m, EigenAlphaObject.getD(), EigenAlphaObject.getV(), nodeRoot.getHeight() - currentSpHeight).getData();  //Math.exp(-alpha * (nodeRoot.getHeight() - currentSpHeight));	
//			WMat[spNr][intRootRegime] += cellValue;
//			
//			while(!sp.isRoot()) {
//					
//				int intRegime = regimes[sp.getNr()].intValue() + rootIndexOffset;	// Will specify the column index
//				cellValue = Math.exp(-alpha * (sp.getHeight() - currentSpHeight)) - Math.exp(-alpha * (sp.getParent().getHeight() - currentSpHeight));
//				WMat[spNr][intRegime] += cellValue;
//					
//				sp = sp.getParent();
//			}
//		}
//		
//		return ;
//	}
	
		
}
