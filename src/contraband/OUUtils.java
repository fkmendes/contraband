package contraband;

import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.lang.Double;


import beast.evolution.tree.Node;
import beast.util.TreeParser;

public class OUUtils {
	
	/*
	 * The math inside this formula comes from 
	 * "Asymptotic theory with hierarchical autocorrelation: Ornstein-Uhlenbeck tree models"
	 * by Lam Si Tun Ho and Cécile Ané (2013)
	 * 
	 * Equation 1 (covariance matrix, random root formula)
	 * 
	 * \frac{\sigma^2}{2\alpha}\mathb{V}\text{ with }V_{ij} = e^{-\alpha d_{ij}}
	 * 
	 * where d_ij is the length of the path between both species
	 * and sigma^2/(2*alpha) is referred to as gamma and is the variance
	 * of the stationary distribution of the process. So the stationary
	 * distribution depends on sigma^2 and alpha (it follows those values
	 * during MCMC). This equation assumes there is a single mean mu (is this
	 * mu supposed to be root optimum?) over the whole tree.
	 * I do not understand how/why this assumption is made since
	 * we are modelling OU in the first place, and in OU you have multiple
	 * optima.
	 * 
	 * Equation 2 (covariance matrix, fixed root formula)
	 * 
	 * \frac{\sigma^2}{2\alpha}\mathb{V}\text{ with }V_{ij} = e^{-\alpha d_{ij}}(1-e^{-2\alpha t_{ij}})
	 * 
	 * where t_{ij} is the covariance between species i and j from the species tree
	 * 
	 * NOTE 1: sigma^2 does not go in because it is handled by MVNUtils
	 * 
	 * NOTE 2: The OU T matrix depends on the "normal" T matrix that comes straight from the tree.
	 * But when we have OU, now this T matrix from the tree has to be modified, which is what
	 * the code below does.
	 */ 
	public static void computeOUTMatOneTrait(int n, double alpha, double[][] tMat, double[][] ouTMat, boolean rootIsFixed) {
		
		double cellValue;
		double divByTwoAlpha = 1.0 / (2.0 * alpha);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				
				cellValue = divByTwoAlpha * Math.exp(-alpha * (tMat[i][i] + tMat[j][j] - 2.0 * tMat[i][j]));
				
				if (rootIsFixed) {
					cellValue *= (1.0 - Math.exp(-2.0 * alpha * tMat[i][j]));
				}
				
				ouTMat[i][j] = cellValue;  // exponent part
			}
		}	
	}
	
	// authors: Marguerite A. Butler* and Aaron A. King
	// title: Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution (Appendix 1)
		// Equation A7 (weight matrix, isolated root formula)
	
	// authors: Julien Clavel, Gilles Escarguel and Gildas Merceron
	// title: mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data
		// Page 5 (weight matrix, merged root idea)
	
	public static void computeWMatOneTrait(TreeParser Tree, int n, int r, double alpha, double[][] WMat, boolean mergeRoot) {
		
		Node nodeRoot = Tree.getRoot();
		List<Node> allLeaf = nodeRoot.getAllLeafNodes();
		
		Double[] regimes = new Double[nodeRoot.getNr() + 1]; // Will keep the vector with the regimes of the branches subtending each node
		int rootIndexOffset;
		int intRootRegime;	// Eldest regime index

		Tree.getMetaData(nodeRoot, regimes, "Regime");
		
		if (mergeRoot) { 		// column dimension of WMat must be r
			
			rootIndexOffset = 0;
			intRootRegime = regimes[nodeRoot.getNr()].intValue();
			
		} else { 		// column dimension of WMat must be r + 1
			
			rootIndexOffset = 1;
			intRootRegime = 0;	
		}
		
		for (Node sp: allLeaf) {
				
			int spNr = sp.getNr();	// Will specify the row index 
				
			// Adding the root chunk in the column of the eldest regime
			double currentSpHeight = sp.getHeight();
			double cellValue = Math.exp(-alpha * (nodeRoot.getHeight() - currentSpHeight));	
			WMat[spNr][intRootRegime] += cellValue;
			
			while(!sp.isRoot()) {
					
				int intRegime = regimes[sp.getNr()].intValue() + rootIndexOffset;	// Will specify the column index
				cellValue = Math.exp(-alpha * (sp.getHeight() - currentSpHeight)) - Math.exp(-alpha * (sp.getParent().getHeight() - currentSpHeight));
				WMat[spNr][intRegime] += cellValue;
					
				sp = sp.getParent();
			}
		}		
		//GeneralUtils.displayRealMatrix(WMat);
	}
	
	
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
