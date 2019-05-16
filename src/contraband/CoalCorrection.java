package contraband;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.TreeParser;

public class CoalCorrection extends CalculationNode {

	final public Input<TreeParser> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<RealParameter> popSizesInput = new Input<>("popSizes", "Temporary input containing the pop size of each node in the tree.", Validate.REQUIRED);
	
	private TreeParser tree;
	private int nSpp;
	
	private double[][] nLineageDistAtEnd; // note that end here means earlier in "regular" (e.g., forward) time
	private double[][] correctedPhyloTMat;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		nSpp = tree.getLeafNodeCount();
		correctedPhyloTMat = new double[nSpp][nSpp];
	}

	private void fillNLineageDistInPlace(Node aNode, Double[] allPopSizes) {
		int nodeIdx = aNode.getNr();
		Double thisNodePopSize = allPopSizes[nodeIdx];

		if (aNode.isLeaf()) {
			nLineageDistAtEnd[nodeIdx] = new double[] { 1.0 };
		}
		
		// if internal node or root
		else {
			Node leftChild = aNode.getChild(0);
			fillNLineageDistInPlace(leftChild, allPopSizes); // recursion
		
			Node rightChild = aNode.getChild(0);
			fillNLineageDistInPlace(rightChild, allPopSizes); // recursion
			
			int maxNLineages = nLineageDistAtEnd[leftChild.getNr()].length + nLineageDistAtEnd[rightChild.getNr()].length;
			nLineageDistAtEnd[nodeIdx] = new double[maxNLineages];
			
			
		}
	}
	
	private void fillPhyloTMatInPlace() {
		tree = treeInput.get();
		Double[] allPopSizes = popSizesInput.get().getValues();
		
		fillNLineageDistInPlace(tree.getRoot(), allPopSizes);
		
		// TODO: then get expected times
		
		// TODO: then get mean total tree heights
		
		// TODO: finally, populate correctedPhyloTMAT and return it
	}
	
	// getters
	public double[][] getCorrectedPhyloTMat() {
		fillPhyloTMatInPlace();
		
		return correctedPhyloTMat;
	}
}
