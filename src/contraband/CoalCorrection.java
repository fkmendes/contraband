package contraband;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import cern.colt.Arrays;

public class CoalCorrection extends CalculationNode {

	final public Input<TreeParser> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<RealParameter> popSizesInput = new Input<>("popSizes", "Temporary input containing the pop size of each node in the tree.", Validate.REQUIRED);
	
	private Tree tree;
	private int nSpp;
	
	private double[][] nLineageDistAtEnd; // note that end here means earlier in "regular" (e.g., forward) time
	private double[][] correctedPhyloTMat;
	
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		nSpp = tree.getLeafNodeCount();
		int nNodes = tree.getNodeCount();
		nLineageDistAtEnd = new double[nNodes][];
		correctedPhyloTMat = new double[nSpp][nSpp];
	}

	private int fillNLineageDistInPlace(Node aNode, Double[] allPopSizes) {
		int maxNLineages = 0;
		int nodeIdx = aNode.getNr();
		Double thisNodePopSize = allPopSizes[nodeIdx];

		if (aNode.isLeaf()) {
			nLineageDistAtEnd[nodeIdx] = new double[] { 1.0 }; // always start with 1 lineage per species
			maxNLineages = 1;
		}
		
		// if internal node or root
		else {
			Node leftChild = aNode.getChild(0);
			int leftIdx = leftChild.getNr();
			int maxNLineagesLeft = fillNLineageDistInPlace(leftChild, allPopSizes); // recursion
		
			Node rightChild = aNode.getChild(1);
			int rightIdx = rightChild.getNr();
			int maxNLineagesRight = fillNLineageDistInPlace(rightChild, allPopSizes); // recursion
			
			maxNLineages = maxNLineagesLeft + maxNLineagesRight;
			nLineageDistAtEnd[nodeIdx] = new double[maxNLineages];
			
			for (int nLineages=2; nLineages <= maxNLineages; ++nLineages) {
				double probOfNLineages = 0.0; // going backward in time, we start at nLineages with this prob
				
				for (int nLeft=1; nLeft <= maxNLineagesLeft; ++nLeft) {
					int nRight = nLineages - nLeft;
					
					if (nRight > 0 && nRight <= maxNLineagesRight) {
						// we're iterating over all combinations of nLeft and nRight that comprise values from 2 to maxNLineages
						probOfNLineages += nLineageDistAtEnd[leftIdx][nLeft-1] * nLineageDistAtEnd[rightIdx][nRight-1]; // -1 is offset
					}
				}
				
				// at the root, we stop at just the probOfNLineages at the start, no Tavare coeffient
				if (aNode.isRoot()) {
					nLineageDistAtEnd[nodeIdx][nLineages-1] = probOfNLineages;
				}
				// at other internal nodes, we have nLineages with probOfNLineages at the start, and then
				// these nLineages might coalesce to kLineages according to Tavare's coefficients
				else {
					for (int kLineages=1; kLineages <= nLineages; ++kLineages) {
						nLineageDistAtEnd[nodeIdx][kLineages-1] = probOfNLineages * 
								CoalUtils.getHeledGij(nLineages, kLineages, aNode.getLength(), thisNodePopSize);
					}
				}
			}
			
		}
		
		System.out.println(Arrays.toString(nLineageDistAtEnd[nodeIdx]));
		return maxNLineages; 
	}
	
	private double getExpGenealogyHeight(Double[] allPopSizes) {
		Node root = tree.getRoot();
		int rootIdx = root.getNr();
		
		double expGenealHeight = 0.0;
		double maxNLineages = fillNLineageDistInPlace(root, allPopSizes);
		
		for (int kLineages=1; kLineages <= maxNLineages; ++kLineages) {
			double probOfKLineages = nLineageDistAtEnd[rootIdx][kLineages-1];
			expGenealHeight += probOfKLineages * CoalUtils.getMeanRootHeight(kLineages, allPopSizes[rootIdx]);
		}
		
		return expGenealHeight + root.getHeight();
	}
	
	private void fillPhyloTMatInPlace() {
		if (treeInput.isDirty()) {
			tree = treeInput.get();
		}
		
		Double[] popSizes = popSizesInput.get().getValues();
		
		double expGenealHeight = getExpGenealogyHeight(popSizes);
		System.out.println("Total genealogical height=" + expGenealHeight);
		
		// TODO: then get expected times
		
		// TODO: finally, populate correctedPhyloTMAT and return it
	}
	
	// getters
	public void getCorrectedPhyloTMat() {
		fillPhyloTMatInPlace();
		
//		for (double[] nodeDist: nLineageDistAtEnd) {
//			System.out.println(Arrays.toString(nodeDist));
//		}
	}
}
