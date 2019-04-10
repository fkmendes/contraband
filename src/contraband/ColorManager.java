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
	
	private int nSpp, nNodes;
	private double[][] spColorValuesMatrix;
	private double[] nodeWeightedColorValues;
	private List<Node> leftLeaves;
	private List<Node> rightLeaves;
	
	private TreeParser tree;
	private Double[] colorValues;
	private Integer[] colorAssignments;
	
	// stored stuff
	double[][] storedSpColorValuesMatrix;
	
	@Override
	public void initAndValidate() {
		readInput();
		checkDimensions();
		
		nSpp = tree.getLeafNodeCount();
		nNodes = tree.getNodeCount();
		nodeWeightedColorValues = new double[nNodes];
		spColorValuesMatrix = new double[nSpp][nSpp];
		leftLeaves = new ArrayList<Node>();
		rightLeaves = new ArrayList<Node>();
		
		storedSpColorValuesMatrix = new double[nSpp][nSpp];
	}

	private void readInput() {
		tree = treeInput.get();
		colorValues = colorValuesInput.get().getValues();
		colorAssignments = colorAssignmentInput.get().getValues();
	}
	
	private void checkDimensions() {
		int nBranches = tree.getNodeCount();
		if (nBranches != colorAssignments.length) {
			throw new RuntimeException("The number of color (rates, or optima) parameters does not match the number of branches in the tree.");
		}
	}
	
	private void populateColorValues() {
		fillNodeColorValues(tree.getRoot());
	}
	
	private void fillNodeColorValues(Node aNode) {
		int nodeIdx = aNode.getNr();
		// System.out.println("nodeIdx=" + nodeIdx); // for debugging 
		
		if (aNode.isLeaf()) {		 
			spColorValuesMatrix[nodeIdx][nodeIdx] = nodeWeightedColorValues[nodeIdx]; // populating diagonal entries
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
				spColorValuesMatrix[aLeftLeaf.getNr()][aRightLeaf.getNr()] = nodeWeightedColorValues[nodeIdx];
				spColorValuesMatrix[aRightLeaf.getNr()][aLeftLeaf.getNr()] = nodeWeightedColorValues[nodeIdx];
			}	
		}
		
		fillNodeColorValues(left);
		fillNodeColorValues(right);
				
		return;
	}
	
	public double[][] getSpColorValuesMatrix() {
		populateColorValues();
		
		return spColorValuesMatrix;
	}
	
	public double getNodeColorValue(Node aNode) {
		colorValues = colorValuesInput.get().getValues();
		colorAssignments = colorAssignmentInput.get().getValues();
		
		return colorValues[colorAssignments[aNode.getNr()]];
	}
	
	@Override
	protected boolean requiresRecalculation() {
		return true;
	}
	
	@Override
	public void store() {
		for (int ithRow=0; ithRow<nNodes; ++ithRow) {
			System.arraycopy(spColorValuesMatrix[ithRow], 0, storedSpColorValuesMatrix[ithRow], 0, spColorValuesMatrix[ithRow].length);
		}
		
		super.store();
	}
	
	@Override
	public void restore() {
		double[][] matTmp;
		
		matTmp = spColorValuesMatrix;
		spColorValuesMatrix = storedSpColorValuesMatrix;
		storedSpColorValuesMatrix = matTmp;
		
		super.restore();
	}
}
