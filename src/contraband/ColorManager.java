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

	final public Input<Integer> nTraitsInput = new Input<>("nTraits", "Number of traits.", Validate.REQUIRED);
	final public Input<Integer> nColorsInput = new Input<>("nColors", "Maximum number of colors.", Validate.REQUIRED);
	final public Input<TreeParser> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<RealParameter> colorValuesInput = new Input<>("colorValues", "Real values (e.g., rates, optima, whatever colors represent) associated to each color (all values for 1st trait, then all values for 2nd trait, and so on).", Validate.REQUIRED);
	final public Input<IntegerParameter> colorAssignmentInput = new Input<>("colorAssignments", "Integers representing colors, one per branch.", Validate.REQUIRED);
	
	final public Input<Boolean> coalCorrectionInput = new Input<>("coalCorrection", "Whether or not to do coalescent correction.", Validate.REQUIRED);
	final public Input<IntegerParameter> rootEdgeColorAssignmentInput = new Input<>("rootEdgeColorAssignment", "Integer representing color of root edge.");
	final public Input<Double> rootEdgeLengthInput = new Input<>("rootEdgeLength", "root edge length.", 0.0);
	
	private boolean doCoalCorrection;
	
	private int nSpp, nNodes; // why part of state? because requires traversing tree
	private double[][] spColorValuesMat;
	private double[] nodeWeightedColorValues;
	private List<Node> leftLeaves;
	private List<Node> rightLeaves;
	private String[] spNamesInASCIIBeticalOrTaxonSetOrder;
	
	private TreeParser tree;
	private Double[] colorValues;
	private Integer[] colorAssignments;
	
	private Double rootEdgeColorValue, rootEdgeLength;
	private Integer rootEdgeColorAssignment;
	private double rootEdgeVar;
	
	// stored stuff
	Double[] storedColorValues;
	Integer[] storedColorAssignments;
//	double[][] storedSpColorValuesMat;
	
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
		spNamesInASCIIBeticalOrTaxonSetOrder = new String[nSpp];
		nodeWeightedColorValues = new double[nNodes];
		spColorValuesMat = new double[nSpp][nSpp];
		leftLeaves = new ArrayList<Node>();
		rightLeaves = new ArrayList<Node>();
		
		// stored stuff
		storedColorValues = new Double[colorValues.length];
		storedColorAssignments = new Integer[nNodes];
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
		
		tree = treeInput.get();
		nSpp = tree.getLeafNodeCount();
		nNodes = tree.getNodeCount();
		colorValues = colorValuesInput.get().getValues();
		colorAssignments = colorAssignmentInput.get().getValues();
		
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
			throw new RuntimeException("The number of color (rates, or optima) assignments does not match the number of branches in the tree. Every branch must be assigned a color.");
		}
		
		// when I implement multiple traits in MVN likelihood, I will probably have to initialize nTraitsˆ2 * maxNColors values... for now this is the number of values we need for pruning multiple
		// characters with multiple colors
		Integer nTraits = nTraitsInput.get();
		Integer nColors = nColorsInput.get();
		if ((nTraits * nColors) != colorValues.length) {
			throw new RuntimeException("The number of initialized color values does not match (max # of colors * # of traits).");
		}
	}
	
	private void populateColorValueMat() {	
		fillNodeColorValuesOneTrait(tree.getRoot(), spNamesInASCIIBeticalOrTaxonSetOrder);
		
		// tree gets taller due to ILS, so it adds some variance to all cells in VCV matrix
		if (rootEdgeVar != 0.0) {
			for (int i=0; i<spColorValuesMat.length; ++i) {
				for (int j=0; j<spColorValuesMat[i].length; ++j) {
					spColorValuesMat[i][j] += rootEdgeVar;
				}
			}
		}
		
//		System.out.println("Just populated color matrix. Species order in VCV Mat is:" + Arrays.toString(spNamesInASCIIBeticalOrTaxonSetOrder));
	}
	
	private void fillNodeColorValuesOneTrait(Node aNode, String[] spOrderInTMat) {
		int nodeIdx = aNode.getNr(); 
		
		if (aNode.isLeaf()) {		
//			System.out.println("Leaf " + aNode.getID() + ", nodeIdx=" + nodeIdx); // uncomment to see index of all leaves 
			spColorValuesMat[nodeIdx][nodeIdx] = nodeWeightedColorValues[nodeIdx]; // populating diagonal entries
			spOrderInTMat[nodeIdx] = aNode.getID();
			return;
		}
		
		if (aNode.isRoot()) { nodeWeightedColorValues[nodeIdx] = 0.0; }
		
		// uncomment to see index of internal nodes and which internal node is which
//		System.out.println("Internal node, nodeIdx=" + nodeIdx); // for debugging
//		List<Node> leafNodes = aNode.getAllLeafNodes();
//		for (Node node: leafNodes) {
//			System.out.println("Daughter leaf: " + node.getID());
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
		
		fillNodeColorValuesOneTrait(left, spOrderInTMat);
		fillNodeColorValuesOneTrait(right, spOrderInTMat);
				
		return;
	}
	
	/*
	 * Getters
	 */
	public Integer getNTraits() {
		return nTraitsInput.get();
	}
	
	public Integer getNColors() {
		return nColorsInput.get();
	}
	
	public int getNNodes() {
		return nNodes;
	}
	
	public Double[] getColorValues() {
		// readInputBeforeGetters();
		colorValues = colorValuesInput.get().getValues();
		
		return colorValues;
	}
	
	public Integer[] getColorAssignments() {
		return colorAssignments;
	}
	
	public double[][] getSpColorValuesMatOneTrait() {
		readInputBeforeGetters();
		populateColorValueMat();
		
		return spColorValuesMat;
	}
	
	public double getNodeColorValue(Node aNode, int traitIdx) {
		readInputBeforeGetters();
		
		Integer nTraits = nTraitsInput.get();
		return colorValues[(nTraits * traitIdx) + colorAssignments[aNode.getNr()]];
	}
	
	public String[] getSpNamesInVCVMatOrder() {
		return spNamesInASCIIBeticalOrTaxonSetOrder;
	}
	
	@Override
	protected boolean requiresRecalculation() {
		return true;
	}
	
	@Override
	public void store() {
		System.arraycopy(colorValues, 0, storedColorValues, 0, colorValues.length);
		System.arraycopy(colorAssignments, 0, storedColorAssignments, 0, colorAssignments.length);
			
		super.store();
	}
	
	@Override
	public void restore() {
		Integer[] vecTmpInt;
		Double[] vecTmpDouble;
		
		vecTmpDouble = colorValues;
		colorValues = storedColorValues;
		storedColorValues = vecTmpDouble;
		
		vecTmpInt = colorAssignments;
		colorAssignments = storedColorAssignments;
		storedColorAssignments = vecTmpInt;
		
		super.restore();
	}
}
