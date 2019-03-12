package test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import beast.evolution.tree.Node;
import beast.util.TreeParser;
import contraband.GeneralUtils;
import contraband.MVNUtils;
import contraband.OUUtils;

public class OUmultiTraitGeneralTest {

	public static void main(String[] args) {
		
		// Setting up the tree
		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// Creating phylogenetic covariance matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();

		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
		
		double[][] Alpha = {{1, 0.2},
							{0.2, 1}
							};
		
		int m = Alpha.length;
		
		// exp(-A * t) test
		RealMatrix alphaMat = new  Array2DRowRealMatrix(Alpha);
		EigenDecomposition eigenAlpha = new EigenDecomposition(alphaMat);
		RealMatrix eigenValAlphaDiag = eigenAlpha.getD();
		System.out.println("eigen values Alpha DiagonalMat");
		GeneralUtils.displayRealMatrix(eigenValAlphaDiag);

		System.out.println("caca");
		RealMatrix eigenVecAlpha = eigenAlpha.getV();
		GeneralUtils.displayRealMatrix(eigenVecAlpha);
		System.out.println("cacaT");
		RealMatrix eigenTVecAlpha = eigenAlpha.getVT();
		GeneralUtils.displayRealMatrix(eigenTVecAlpha);
		GeneralUtils.displayRealMatrix(eigenValAlphaDiag);
		RealMatrix expAlphaMat = OUUtils.getExpAlphaMat(m, eigenValAlphaDiag, eigenVecAlpha, 1);
		GeneralUtils.displayRealMatrix(eigenValAlphaDiag);
		
		System.out.println("Alpha diagonal Mat");
		GeneralUtils.displayRealMatrix(eigenValAlphaDiag);
		System.out.println("exp trial");
		GeneralUtils.displayRealMatrix(expAlphaMat);
		
		
		// computeExpAlphaMat works for symmetric Alpha matrices
		
		// mainChunk test
		System.out.println();
		GeneralUtils.display2DArray(tMatInput);
		double[][] mainChunk = new double[m][m];
		OUUtils.treeAlphaChunkMat(m, eigenAlpha.getRealEigenvalues(), tMatInput[0][0], tMatInput[1][1], tMatInput[0][1], true, mainChunk);
		System.out.println("mainChunk");
		GeneralUtils.display2DArray(mainChunk);
		// mainChunkMat works fine
		
		//eigencholChunk test
		double[][] Sigma2 = {{2, 0.5},
							{0.5, 3}
							};
		RealMatrix Sigma2Mat = new Array2DRowRealMatrix(Sigma2);
		CholeskyDecomposition cholSigma2 = new CholeskyDecomposition(Sigma2Mat);
		GeneralUtils.displayRealMatrix(cholSigma2.getL().multiply(cholSigma2.getLT()));
		
		RealMatrix ecChunk = OUUtils.eigencholChunk(cholSigma2.getL(), eigenVecAlpha);
		System.out.println("ecChunk");
		GeneralUtils.displayRealMatrix(ecChunk);
		
		// I think it works well despite symmetry
		
		// Element wise product
		double[][] elementWiseRes = new double[m][m];
		GeneralUtils.elementWiseProduct(mainChunk, ecChunk.getData(), elementWiseRes);
		System.out.println("elementwise");
		GeneralUtils.display2DArray(elementWiseRes);
		
		RealMatrix elementWiseMat = new Array2DRowRealMatrix(elementWiseRes);
		System.out.println("bigChunk");
		GeneralUtils.displayRealMatrix(eigenVecAlpha.multiply(elementWiseMat).multiply(eigenTVecAlpha));
		
		// Results match with R

		// Let's see the covariance matrix between species 1 and 2
		RealMatrix cov12 = OUUtils.multiTraitCovMat2Species(m, ecChunk.getData(), eigenAlpha, 3, 3, 2, true);
		System.out.println("This is covariance of sp 1 and 2, fixedRoot");
		GeneralUtils.displayRealMatrix(cov12);
		// Works well for the fixedRoot case
		
		// See the randomRoot case
		mainChunk = new double[m][m];
		OUUtils.treeAlphaChunkMat(m, eigenAlpha.getRealEigenvalues(), tMatInput[0][0], tMatInput[1][1], tMatInput[0][1], false, mainChunk);
		GeneralUtils.display2DArray(mainChunk);
		cov12 = OUUtils.multiTraitCovMat2Species(m, ecChunk.getData(), eigenAlpha, 3, 3, 2, false);
		System.out.println("This is covariance of sp 1 and 2, randomRoot");
		GeneralUtils.displayRealMatrix(cov12);
		// Works fine as well
		
		// weight Matrix sp i Vs regime k
		Node nodeRoot = myTree.getRoot();
		List<Node> allLeaf = nodeRoot.getAllLeafNodes();
		Double[] regimes = new Double[nodeRoot.getNr() + 1]; //It stores the regimes of nodes in the order taken by beast
		myTree.getMetaData(nodeRoot, regimes, "Regime");	// Populating regimes vector
		Node sp = allLeaf.get(0);
		System.out.println("species index: " + sp.getNr());
		
//		int intRootRegime = regimes[nodeRoot.getNr()].intValue(); //We will set the regime of the root to be always 0
		
		RealMatrix W01 = OUUtils.multiTraitWeightMatSpVsRegime(m, sp, 1, nodeRoot, regimes, eigenAlpha, true);
		System.out.println("tadaaa");
		GeneralUtils.displayRealMatrix(W01);
		
		// Test a forLoopRange
		int[] Lista = GeneralUtils.forLoopRange(20, 10); //Lista de 20 hasta 29
		System.out.println("Lista de prueba = " + Lista[2]);
		
		// Works fine
		
		// Compute the full weight matrix
		int r = 3; //number of regimes
		double[][] Result = new double [n * m][(r + 1) * m];
		OUUtils.multiTraitFullWeightMat(m, nodeRoot, regimes, 3, eigenAlpha, false, Result);
		GeneralUtils.display2DArray(Result);
		
	}
}
