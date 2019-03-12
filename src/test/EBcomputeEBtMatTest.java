package test;
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import beast.evolution.tree.Node;
import beast.util.TreeParser;
import contraband.EBUtils;
import contraband.GeneralUtils;
import contraband.MVNUtils;


// In this test we will verify the value of function GChunk, kronecker and computeEBtMat for an Early Burst model.

// 1, 2, 3 refers to Test 1 and Test 2 respectively (WE ARE USING THE FIRST TWO TREES AS IN OUOneTraitcomputeOUTMatOneTraitTest.java)
// We will also show how to compute the likelihood in this model

public class EBcomputeEBtMatTest {
	
	final static double EPSILON = 1e-4;
	private static double[][] res1, res2;
	
	@BeforeClass	// Test 1: Ultrametric tree
	public static void setUP1() throws Exception {
		
		String treeStr = "(((sp1:1.0, sp2:1.0):1.0, sp3:2.0):1.0, sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// Data input for covariance function
		double[][] sigmaMat = { {0.003311920, 0.007755767},		// hardcoded from r mvEB result
				{0.007755767, 0.033005583}
				};
		int m = sigmaMat.length;
		//double[] theta = {-0.2417452, -0.3869768};
		double g = 1.999999;
		
		// Creating phylogenetic covariance matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();

		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
		
		// GChunk and covariance calculation
		double[][] gMat = new double[n][n];
		res1 = new double[n * m][n * m];
		EBUtils.Gchunk(g, tMatInput, gMat);
		
		EBUtils.computeEBtMat(gMat, sigmaMat, res1);
		
	}
	
	@BeforeClass	// Test 2: nonUltrametric tree
	public static void setUP2() throws Exception {
		
		String treeStr = "(((sp1:2.0, sp2:1.0):1.0, sp3:4.0):1.0, sp4:3.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// Data input for covariance function
		double[][] sigmaMat = { {0.4782103, -0.6159709},		// hardcoded from r mvEB result
								{-0.6159709, 0.8544459}
								};
		int m = sigmaMat.length;
		//double[] theta = {0.5941884, -1.261628};
		double g = 0.2204587;
		
		// Creating phylogenetic covariance matrix of the tree
		int n = myTree.getLeafNodeCount();
		double[][] tMatInput = new double[n][n];
		double[] nodeToRootPaths = new double[myTree.getNodeCount()];
		List<Node> leftLeaves = new ArrayList<>();
		List<Node> rightLeaves = new ArrayList<>();

		MVNUtils.populateTMatrix(myTree, nodeToRootPaths, tMatInput, leftLeaves, rightLeaves);
		
		// GChunk and covariance calculation
		double[][] gMat = new double[n][n];
		res2 = new double[n * m][n * m];
		EBUtils.Gchunk(g, tMatInput, gMat);
		
		EBUtils.computeEBtMat(gMat, sigmaMat, res2);
		
	}
	
	
	
	@Test // Test 1
	public void againstRcovMatEBTest1 () {
		
		Assert.assertEquals(0.66640428, res1[0][0], EPSILON);
		Assert.assertEquals(0.08875625, res1[0][2], EPSILON);
		Assert.assertEquals(0.2078471, 	res1[2][1], EPSILON);
	}
	
	@Test // Test 2
	public void againstRcovMatEBTest2 () {
		
		Assert.assertEquals(3.0700716, res2[0][0], EPSILON);
		Assert.assertEquals(1.2020018, res2[0][2], EPSILON);
		Assert.assertEquals(-1.5482690,res2[2][1], EPSILON);	
	}

}
