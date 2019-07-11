package contraband;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeUtils;
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

	/*
	 * Does not assume tree is ultrametric
	 */
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
	
	/*
	 *  Note that this is the expected genealogy height where the "tips" start at the root population
	 *  so this is NOT the TOTAL genealogy height (as in, the total gene tree height)
	 */
	private double getExpGenealHeightAtRoot(Double[] allPopSizes) {
		Node root = tree.getRoot();
		int rootIdx = root.getNr();
		
		double expGenealHeight = 0.0;
		double maxNLineages = fillNLineageDistInPlace(root, allPopSizes);
		
		for (int kLineages=1; kLineages <= maxNLineages; ++kLineages) {
			double probOfKLineages = nLineageDistAtEnd[rootIdx][kLineages-1];
			expGenealHeight += probOfKLineages * CoalUtils.getMeanRootHeight(kLineages, allPopSizes[rootIdx]);
		}
		
		return expGenealHeight;
	}
	
	/*
	 * This works for non-ultrametric trees as well, the only thing that matters
	 * is the height of the MRCA of the two pairs of species, and expected times
	 * are all relative to the most recent leaf in the tree, which is what we end
	 * up using in the calculation of variances later
	 */
	private double getExpCoalTimePair(Node node1, Node node2, Double[] allPopSizes) {
		Set<String> pairOfSpNames = new HashSet<String>();
		pairOfSpNames.add(node1.getID());
		pairOfSpNames.add(node2.getID());
		
		Node mrcaInSpTree = TreeUtils.getCommonAncestorNode(tree, pairOfSpNames); // mrca of two nodes in species tree
		double mrcaHeightSpTree = mrcaInSpTree.getHeight();
		
		double expCoalTimePair = 0.0;
		double probCurr = 1.0;
		
		while(true) {
			double branchLength = mrcaInSpTree.getLength();
			double popSize = allPopSizes[mrcaInSpTree.getNr()];
			double invPopSize = 1.0 / popSize;
			
			double e = 0.0;
			double probNoCoal = 0.0;
			
			if (mrcaInSpTree.isRoot()) {
				e = popSize;
				probNoCoal = 0.0;
			}
			
			else {
				e = (1 - (branchLength * invPopSize + 1) * Math.exp(-branchLength * invPopSize)) / invPopSize;
				probNoCoal = CoalUtils.getHeledGij(2, 2, branchLength, popSize);
			}
			
			// expCoalTimePair += (mrcaHeightSpTree + e) * probCurr * (1-probNoCoal); // original was bugged
			expCoalTimePair += (mrcaHeightSpTree + e/(1-probNoCoal)) * (1-probNoCoal) * probCurr;
			
			if (mrcaInSpTree.isRoot()) {
				break;
			}
			
			probCurr *= probNoCoal;
			mrcaHeightSpTree += branchLength;
			mrcaInSpTree = mrcaInSpTree.getParent();
		}
		
		return expCoalTimePair;
	}
	
	private void fillPhyloTMatInPlace() {
		if (treeInput.isDirty()) {
			tree = treeInput.get();
		}
		
		double treeHeight = tree.getRoot().getHeight();
		List<Node> allLeaves = tree.getRoot().getAllLeafNodes();
		Double[] popSizes = popSizesInput.get().getValues();
		
		double expGenealHeightAtRoot = getExpGenealHeightAtRoot(popSizes);
		
		System.out.println("Expected total genealogical height = " + (treeHeight+expGenealHeightAtRoot));
		
		// TODO: then get expected times
		List<Node> allNodes = tree.getRoot().getAllLeafNodes();
		for (int i=0; i<nSpp; ++i) {
			for (int j=0; j<nSpp; ++j) {
				if (i == j) {
					correctedPhyloTMat[i][j] = expGenealHeightAtRoot + (treeHeight - allLeaves.get(i).getHeight()); // all variances (diag elements) don't have to be the same
				}
				
				else {
					// we subtract the expected coalescent times from the total genealogical height (= root height + expGenealHeightAtRoot)
					correctedPhyloTMat[i][j] = tree.getRoot().getHeight() + expGenealHeightAtRoot - getExpCoalTimePair(allNodes.get(i), allNodes.get(j), popSizes);
					correctedPhyloTMat[j][i] = correctedPhyloTMat[i][j]; // assume matrix is symmetric
				}
			}
		}
	}
	
	// getters
	public double[][] getCorrectedPhyloTMat() {
		fillPhyloTMatInPlace();
		
		GeneralUtils.display2DArray(correctedPhyloTMat);
		
		return correctedPhyloTMat;
	}
}
