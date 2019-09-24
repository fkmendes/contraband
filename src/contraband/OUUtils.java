package contraband;

import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.lang.Double;

import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import cern.colt.Arrays;

/*
 * Notes about Hansen model (=OU on a tree with multiple regimes)
 * 
 * In OU, there will be r selective optima or regimes. theta_k is the k-th regime. 
 * The total number of regimes, r, is independent of the number of branches in the tree
 * (and will usually be much much smaller), but theta_0 is the separate regime (optimum)
 * the root might or might not be assigned.
 * 
 * Note that the likelihood described in Butler & King is GIVEN y_0, the root value
 * (which is often set to be = to the root optimum, theta_0 -- hence why the root state
 * is often referred to as theta_0). As C. Ané pointed out to me, this is a very liberal assumption.
 * In other words, the process is conditioned on a particular y_0 (theta_0) that we can sample during MCMC,
 * for example). In other words, we sample theta_0, then fix it for the rest of the process.
 *
 * However, we can also assume y_0 has a stationary distribution with mean theta_0 and
 * variance = \gamma = \frac{\sigma}{2\alpha}. Think of this stationary distribution as a parametric (depends on the actual parameters
 * you are trying to estimate) prior.
 * This changes the variance of the OU MVN likelihood (what here we call the OU T matrix), as we
 * need to add a bit more variance due to treating theta_0 as a random variable.
 * The mean of the OU process is still computed the same way, as it depends on what theta0 is.
 * 
 * The expected values at any time point in the tree (ancestral states) will be given by the
 * weight matrix W (i.e., the design matrix, which depends on alpha) * the theta vector (or you
 * can just compute the spelled-out formula in Butler & King + others). Doing W * theta vector is a trick
 * that uses the GLS formulation. Differently from BM, the root value is not the expected value for each
 * lineage (the latter being W * theta vector).
 * Note that the thetas are not the ancestral states, they are just the optima pulling the trait values.
 * In the specific case of the root, the root optimum (theta_0) will be the root value (i.e., the mean of
 * the process = vector of theta_0s) if we assume stationarity (rootIsFixed=false). 
 * 
 * ---
 * 
 * Notes about this implementation
 * 
 * While technically you can do mergeRoot=false/true with rootIsFixed=false/true in any
 * combination you like, it makes sense to only do mergeRoot=false with rootIsFixed=true,
 * and only with non-ultrametric trees (with fossils). If you can't know anything about
 * the root state you probably can't estimate its optimum; so rootIsFixed=false with
 * mergeRoot=true probably isn't worth the additional parameter time. In some packages,
 * the default is actually mergeRoot=true as most trees are ultrametric.
 */

public class OUUtils {
	
	/*
	 * The math inside this formula comes from 
	 * "Asymptotic theory with hierarchical autocorrelation: Ornstein-Uhlenbeck tree models"
	 * by Lam Si Tun Ho and Cécile Ané (2013)
	 * 
	 * Equation 1 (OU T matrix, random root formula)
	 * 
	 * \frac{\sigma^2}{2\alpha}\mathb{V}\text{ with }V_{ij} = e^{-\alpha d_{ij}}
	 * 
	 * where d_ij is the length of the path between both species
	 * and sigma^2/(2*alpha) is referred to as gamma and is the variance
	 * of the Gaussian stationary distribution of the process.
	 * 
	 * This is for ultrametric trees, where there isn't information about the root (mean) values.
	 * This formula (="rootIsFixed=false") assumes we integrate over a stationary distribution with mean=theta_0 (i.e., the
	 * first element of the theta vector) and variance is sigma/(2*alpha).
	 * 
	 * Equation 2 (OU T matrix, fixed root formula); this also appears in the appendix of Butler and King as 
	 * equation A5
	 * 
	 * \frac{\sigma^2}{2\alpha}\mathb{V}\text{ with }V_{ij} = e^{-\alpha d_{ij}}(1-e^{-2\alpha t_{ij}})
	 * 
	 * where t_{ij} is the covariance between species i and j from the species tree.
	 * 
	 * This is for non-ultrametric trees, where fossils are available -- and the root state
	 * is not assumed to come from a distribution centered at theta0's, but are conditioned upon.
	 * But if the tree is ultrametric and OU is stationary, this option should converge on the random root case.
	 * 
	 * NOTE 1: sigma^2 does not go in because it is handled by MVNUtils.
	 * 
	 * NOTE 2: The OU T matrix depends on the "normal" T matrix that comes straight from the tree.
	 * But when we have OU, now this T matrix from the tree has to be modified, which is what
	 * the code below does. It is not the vcv matrix yet because we haven't multiplied it
	 * by sigma^2.
	 */ 
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
	
	/*
	 * The math inside this formula comes from 
	 * 
	 * "Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution" (Appendix 1)
	 * by Marguerite Butler and Aaron King (2004)
	 * 
	 * Equation A7 (weight matrix, root has its own optimum, i.e., mergeRoot=false)
	 * 
	 * W_{ij} = e^{-\alpha T}\sum_{\gamma =1}^{\kappa (i)}\beta_{ik}^{\gamma}(e^{\alpha t_{i}^{\gamma}}-e^{\alpha t_{i}^{\gamma -1}})
	 * 
	 * where T is the length of the time the regime applies
	 * \beta_{ik} is the indicator value of the k-th regime (out of r), i-th species (out of N) (there is one \beta_{ik} per branch)
	 * \gamma is a discrete variable that represents branches on a path (\gamma - 1 is the parent of \gamma), with \gammas going from
	 * 1 to the total number of branches on a path, \kappa{i}.
	 * 
	 * This is the more general description of the weight matrix that appears in Butler & King, in which we allow the root to be annotated
	 * in the tree and have its own optimum (theta_0).
	 * We should use this when we have fossils (tree is non-ultrametric), and this optimum can be estimated.
	 * 
	 * and
	 * 
	 * "mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data"
	 * by Julien Clavel, Gilles Escarguel and Gildas Merceron (2015)
	 * 
	 * Equation 6 (they call it the covariance matrix, this is the multivariate equivalent of the weight matrix in Butler and King)
	 * 
	 * \text{Cov}(Y_i,Y_j)=\int_{0}^{C_{ij}}e^{-\mathbf{A}(C_{ii}-\nu)}\mathbf{R}e^{-\mathbf{A}^{T}(C_{jj}-\nu)}\text{d}\nu
	 * 
	 * In this less general version, we ignore the root optimum metadata, and assume theta_0 is the same as one of the other thetas
	 * (one of the optima). We should use this with ultrametric trees.
	 */
	 public static void computeWMatOneTrait(Integer[] allNodeRegimes, Node rootNode, List<Node> allLeafNodes, int n, int r, double alpha, RealMatrix wMat, boolean useRootMetaData) {		
		 int rootIndexOffset;
		 int rootRegimeIdx;	// Specified in colorAssignment (last position)
			
		 if (useRootMetaData) {
			 rootIndexOffset = 1; // column dimension of WMat must be r + 1
			 rootRegimeIdx = allNodeRegimes[rootNode.getNr()].intValue(); // getting meta data!
		 }
		 else {
			 rootIndexOffset = 0; // column dimension of WMat must be r	
			 rootRegimeIdx = 0; // not getting metadata, and using the color that has index 0
		 }
		 
		 for (Node sp: allLeafNodes) {	
			 int spNr = sp.getNr();	// Will specify the row index 

			 // Adding the root chunk in the column of the eldest regime
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
	 public static void populateOUMeanVector(double alpha, double rootValue, Node aNode, Node rootNode, List<Node> allLeafNodes, BranchRateModel clock, RealVector ouMeanVector, boolean useRootMetaData, double avgSoFar) {
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
//		 double runningAvg = avgSoFar;
//
//		 if (!aNode.isRoot()) {
//			 System.out.println("parent height=" + aNode.getParent().getHeight());
//			 System.out.println("child height=" + aNode.getHeight());
//			 runningAvg += (Math.exp(-alpha * aNode.getParent().getHeight()) - Math.exp(-alpha * aNode.getHeight())) * clock.getRateForBranch(aNode);
//		 }
//
//		 for (Node descNode : aNode.getChildren()) {
//			 populateOUMeanVector(alpha, rootValue, descNode, rootNode, allLeafNodes, clock, ouMeanVector, useRootMetaData, runningAvg);
//		 }
//
//		 if (aNode.isLeaf()) {
//			if (useRootMetaData) {
//				runningAvg += (Math.exp(-alpha * (rootNode.getHeight() - aNode.getHeight()))) * rootValue;
//			}
//
//			ouMeanVector.setEntry(aNode.getNr(), runningAvg);
//		}

		 /*
		  * This implementation below follows Eq. 3 in Hansen 1997.
		  * 
		  * But note that it traverses the tree n times, where
		  * n is the number of species. It is slower than the
		  * recursive solution above, but works for both ultrametric
		  * and non-ultrametric trees
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
				 cellValue = Math.exp(-alpha * tAi) * clock.getRateForBranch(rootNode); // use theta_0
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
	 
	 /*
	  * DEPRECATED: uses tree meta data, cannot operate on it
	  */
//	 public static void computeWMatOneTrait2(TreeParser tree, Node rootNode, List<Node> allLeafNodes, int n, int r, double alpha, RealMatrix wMat, boolean useRootMetaData) {		
//			Double[] allNodeRegimes = new Double[rootNode.getNr() + 1]; // Will keep the vector with the regimes of the branches subtending each node
//			int rootIndexOffset;
//			int rootRegimeIdx;	// Eldest regime index
//
//			/* 
//			 * The regimes array is the computational equivalent of 
//			 * Beta_{i}^{\gamma} in equation A6 of Butler & King's Appendix 1
//			*/
//			tree.getMetaData(rootNode, allNodeRegimes, "Regime"); // writes on regimes array
//			
//			if (!useRootMetaData) { 		// column dimension of WMat must be r			
//				rootIndexOffset = 0;
//				rootRegimeIdx = allNodeRegimes[rootNode.getNr()].intValue();
//				
//			} else { 		// column dimension of WMat must be r + 1
//				rootIndexOffset = 1;
//				rootRegimeIdx = 0;	
//			}
//			
//			for (Node sp: allLeafNodes) {	
//				int spNr = sp.getNr();	// Will specify the row index 
//					
//				// Adding the root chunk in the column of the eldest regime
//				double currentSpHeight = sp.getHeight();
//				double cellValue = Math.exp(-alpha * (rootNode.getHeight() - currentSpHeight));	
//				wMat.addToEntry(spNr, rootRegimeIdx, cellValue);
//				
//				while(!sp.isRoot()) {
//					int regimeIdx = allNodeRegimes[sp.getNr()].intValue() + rootIndexOffset;	// Will specify the column index
//					cellValue = Math.exp(-alpha * (sp.getHeight() - currentSpHeight)) - Math.exp(-alpha * (sp.getParent().getHeight() - currentSpHeight));
//					wMat.addToEntry(spNr, regimeIdx, cellValue);
//						
//					sp = sp.getParent();
//				}
//			}		
//		}
	 
	// MultiTrait Covariance matrix functions
	
	public static RealMatrix getExpAlphaMat(int m, RealMatrix eigenDiag, RealMatrix eigenVecMat, double time) {
		
		RealMatrix Result = eigenDiag.createMatrix(m, m);
		
		for (int i = 0; i < m; i++) {
			Result.setEntry(i, i, Math.exp( - time * eigenDiag.getEntry(i, i)));
		}
		
		return ( Result
				.preMultiply(eigenVecMat)
				.multiply(eigenVecMat.transpose()) );
	}
	
	public static void treeAlphaChunkMat(int m, double[] eigenValues, double cii, double cjj, double cij, boolean fixedRoot, double[][] Result) {
		
		if(fixedRoot) {
			
			for (int i = 0; i < m; i++) {
				
			    for (int j = 0; j < m; j++) {	
			    	
			    	Result[i][j] = ( 1 - Math.exp( - (eigenValues[i] + eigenValues[j]) * cij ) )
			    			/( eigenValues[i] + eigenValues[j] );   
			    } 
			}
			
		} else {
			
			for (int i = 0; i < m; i++) {
				
			    for (int j = 0; j < m; j++) {	
			    	
			    	Result[i][j] = ( Math.exp( - eigenValues[i] * (cii - cij) ) 
			    			* Math.exp( - eigenValues[j] * (cjj -cij) ) )/( eigenValues[i] + eigenValues[j] );   
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
			GeneralUtils.elementWiseProduct(res, eigenCholChunk, res);
			resMat = new Array2DRowRealMatrix(res);

			resMat = OUUtils.getExpAlphaMat(m, eigenObjectAlpha.getD(), eigenObjectAlpha.getV(), cii - cij)
					.multiply(eigenObjectAlpha.getV())
					.multiply(resMat)
					.multiply(eigenObjectAlpha.getVT())
					.multiply(OUUtils.getExpAlphaMat(m, eigenObjectAlpha.getD(), eigenObjectAlpha.getV(), cjj - cij));
		} else {
			
			OUUtils.treeAlphaChunkMat(m, eigenObjectAlpha.getRealEigenvalues(), cii, cjj, cij, false, res);
			GeneralUtils.elementWiseProduct(res, eigenCholChunk, res);
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
					GeneralUtils.setBlocksInMatrix(blockWeightMat.getData(), m, sp.getNr(), k, Result);
				}
			}
			
		} else {
			
			for(Node sp:allLeaf) {
				
				double spHeight = sp.getHeight();
				RealMatrix cellValue = OUUtils.getExpAlphaMat(m, EigenAlphaObject.getD(), EigenAlphaObject.getV(), nodeRoot.getHeight() - spHeight); 
				GeneralUtils.setBlocksInMatrix(cellValue.getData(), m, sp.getNr(), 0, Result);	
				
				for(int k = 0; k < r; k++) {
		
					blockWeightMat = multiTraitWeightMatSpVsRegime(m, sp, k, nodeRoot, regimes, EigenAlphaObject, false);
					GeneralUtils.setBlocksInMatrix(blockWeightMat.getData(), m, sp.getNr(), k + 1, Result);	
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
