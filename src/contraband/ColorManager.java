package contraband;

import java.util.ArrayList;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.TreeParser;

/*
 * This class implements an uncorrelated discrete-distribution local color
 * where color here can be the evolutionary rate, or OU optimum 
 */
public class ColorManager extends CalculationNode {

	final public Input<TreeParser> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<RealParameter> colorValuesInput = new Input<>("colorValues", "Real values (e.g., rates, optima, whatever colors represent) associated to each color.", Validate.REQUIRED);
	final public Input<IntegerParameter> colorAssignmentInput = new Input<>("colorAssignments", "Integers representing colors, one per branch.", Validate.REQUIRED);
	
	final public Input<Boolean> coalCorrectionInput = new Input<>("coalCorrection", "Whether or not to do coalescent correction.", Validate.REQUIRED);
	final public Input<IntegerParameter> rootEdgeColorAssignmentInput = new Input<>("rootEdgeColorAssignment", "Integer representing color of root edge.");
	final public Input<Double> rootEdgeLengthInput = new Input<>("rootEdgeLength", "root edge length.", 0.0);
	
	private boolean doCoalCorrection;
	
	private int nSpp, nNodes;
	private double[][] spColorValuesMat;
	private double[] nodeWeightedColorValues;
	private List<Node> leftLeaves;
	private List<Node> rightLeaves;
	private String[] spNamesInPhyloTMatOrder;
	
	private TreeParser tree;
	private Double[] colorValues;
	private Integer[] colorAssignments;
	
	private Double rootEdgeColorValue, rootEdgeLength;
	private Integer rootEdgeColorAssignment;
	private double rootEdgeVar;
	
	// stored stuff
	double[][] storedSpColorValuesMat;
	
	@Override
	public void initAndValidate() {	
		tree = treeInput.get();
		colorValues = colorValuesInput.get().getValues();
		colorAssignments = colorAssignmentInput.get().getValues();

		checkDimensions();
		
		doCoalCorrection = coalCorrectionInput.get();
		
		// TODO: this will change when I link it to the MSC
		if (doCoalCorrection) {
			rootEdgeColorAssignment = rootEdgeColorAssignmentInput.get().getValue();
			rootEdgeColorValue = colorValues[rootEdgeColorAssignment];
			rootEdgeLength = rootEdgeLengthInput.get();
			rootEdgeVar = rootEdgeColorValue * rootEdgeLength;
		}
		
		nSpp = tree.getLeafNodeCount();
		nNodes = tree.getNodeCount();
		spNamesInPhyloTMatOrder = new String[nSpp];
		nodeWeightedColorValues = new double[nNodes];
		spColorValuesMat = new double[nSpp][nSpp];
		leftLeaves = new ArrayList<Node>();
		rightLeaves = new ArrayList<Node>();
		
		storedSpColorValuesMat = new double[nSpp][nSpp];
	}

	private void readInputBeforeGetters() {
		if (treeInput.isDirty()) {
			tree = treeInput.get();
			
			// TODO: this will change when I link it to the MSC
			if (doCoalCorrection) {
				rootEdgeColorAssignment = rootEdgeColorAssignmentInput.get().getValue();
				rootEdgeColorValue = colorValues[rootEdgeColorAssignment];
				rootEdgeLength = rootEdgeLengthInput.get();
				rootEdgeVar = rootEdgeColorValue * rootEdgeLength;
			}
		}
		
		if (colorValuesInput.isDirty()) {
			colorValues = colorValuesInput.get().getValues();
		}
		
		if (colorAssignmentInput.isDirty()) {
			colorAssignments = colorAssignmentInput.get().getValues();
		}
	}
	
	private void checkDimensions() {
		int nBranches = tree.getNodeCount();
		if (nBranches != colorAssignments.length) {
			throw new RuntimeException("The number of color (rates, or optima) parameters does not match the number of branches in the tree.");
		}
	}
	
	private void populateColorValueMat() {
		fillNodeColorValues(tree.getRoot(), spNamesInPhyloTMatOrder);
		
		// tree gets taller due to ILS, so it adds some variance to all cells in VCV matrix
		if (rootEdgeVar != 0.0) {
			for (int i=0; i<spColorValuesMat.length; ++i) {
				for (int j=0; j<spColorValuesMat[i].length; ++j) {
					spColorValuesMat[i][j] += rootEdgeVar;
				}
			}
		}
	}
	
	private void fillNodeColorValues(Node aNode, String[] spOrderInTMat) {
		int nodeIdx = aNode.getNr();
		// System.out.println("nodeIdx=" + nodeIdx); // for debugging 
		
		if (aNode.isLeaf()) {		 
			spColorValuesMat[nodeIdx][nodeIdx] = nodeWeightedColorValues[nodeIdx]; // populating diagonal entries
			spOrderInTMat[nodeIdx] = aNode.getID();
			return;
		}
		
		if (aNode.isRoot()) { nodeWeightedColorValues[nodeIdx] = 0.0; }
		
		// for debugging
//		List<Node> leafNodes = aNode.getAllLeafNodes();
//		for (Node node: leafNodes) {
//			System.out.println("leaf: " + node.getID());
//		}
		
		Node left = aNode.getChild(0);
		int leftIdx = left.getNr();
		double leftColorValue = colorValues[colorAssignments[leftIdx]];
		nodeWeightedColorValues[leftIdx] = nodeWeightedColorValues[nodeIdx] + (left.getLength() * leftColorValue);
			
		Node right = aNode.getChild(1);
		int rightIdx = right.getNr();
		double rightColorValue = colorValues[colorAssignments[rightIdx]];
		nodeWeightedColorValues[rightIdx] = nodeWeightedColorValues[nodeIdx] + (right.getLength() * rightColorValue);
		
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
				spColorValuesMat[aLeftLeaf.getNr()][aRightLeaf.getNr()] = nodeWeightedColorValues[nodeIdx];
				spColorValuesMat[aRightLeaf.getNr()][aLeftLeaf.getNr()] = nodeWeightedColorValues[nodeIdx];
			}	
		}
		
		fillNodeColorValues(left, spOrderInTMat);
		fillNodeColorValues(right, spOrderInTMat);
				
		return;
	}
	
	public double[][] getSpColorValuesMat() {
		readInputBeforeGetters();
		populateColorValueMat();
		
		return spColorValuesMat;
	}
	
	public double getNodeColorValue(Node aNode) {
		readInputBeforeGetters();
		
		return colorValues[colorAssignments[aNode.getNr()]];
	}
	
	public String[] getSpNamesInPhyloTMatOrder() {
		return spNamesInPhyloTMatOrder;
	}
	
	@Override
	protected boolean requiresRecalculation() {
		return true;
	}
	
	@Override
	public void store() {
		for (int ithRow=0; ithRow<nNodes; ++ithRow) {
			System.arraycopy(spColorValuesMat[ithRow], 0, storedSpColorValuesMat[ithRow], 0, spColorValuesMat[ithRow].length);
		}
		
		super.store();
	}
	
	@Override
	public void restore() {
		double[][] matTmp;
		
		matTmp = spColorValuesMat;
		spColorValuesMat = storedSpColorValuesMat;
		storedSpColorValuesMat = matTmp;
		
		super.restore();
	}
}
