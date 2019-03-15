package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.evolution.tree.Node;
import beast.util.TreeParser;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import contraband.GeneralUtils;
import contraband.MVNUtils;

public class BMOneTraitComputeMVNLk {
	
	final static double EPSILON = 1e-4;
	private double resLk; 
	
	@Before
	public void setUP() throws Exception {
		
		// tree
		String treeStr = "(((sp1:1.0, sp2:1.0):1.0, sp3:2.0):1.0, sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();
			
		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
		GeneralUtils.display2DArray(tMatInput);

		// var, mean
		double var = 1.4822794118;
		double[] meanInput = new double[] { 3.079142, 3.079142, 3.079142, 3.079142 };
		RealVector mean = new ArrayRealVector(meanInput);
				
		// trait values
		double[] dataInput = new double[] { 4.1, 4.5, 5.9, 0.0 };
		RealVector data = new ArrayRealVector(dataInput);
					
		RealMatrix tMat = new Array2DRowRealMatrix(tMatInput);
		LUDecomposition tMatLUD = new LUDecomposition(tMat);
		RealMatrix invTMat = tMatLUD.getSolver().getInverse(); // if only variance changes and not tree, we don't have to invert
		RealMatrix invVcvMat = invTMat.scalarMultiply(1/var); // (var*tMat)^-1 = tMat^-1 / var
		double varToNdetTMat = Math.pow(var, n) * tMatLUD.getDeterminant(); // det(var*tMat) = var^n * det(tMat)
		
		resLk = MVNUtils.getMVNLk(n, mean, data, invVcvMat, varToNdetTMat);
	}
	
	@Test 
	public void againstresLikelihood () {
		Assert.assertEquals(0.0002498382, resLk, EPSILON);
	}


}
