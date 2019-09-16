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

public class CoalCorrection extends CalculationNode {

	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<RealParameter> popSizesInput = new Input<>("popSizes", "Temporary input containing the pop size of each node in the tree.", Validate.REQUIRED);
	
	private Tree tree;
	
	private double[][] nLineageDistAtEnd; // note that end here means earlier in "regular" (e.g., forward) time
	private double[][] correctedPhyloTMat;

	// stored stuff
	private double[][] storedCorrectedPhyloTMat;
	private double[][] storedNLineageDistAtEnd;

	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		int nSpp = tree.getLeafNodeCount();
		int nNodes = tree.getNodeCount();
		nLineageDistAtEnd = new double[nNodes][];
		correctedPhyloTMat = new double[nSpp][nSpp];

		storedNLineageDistAtEnd = new double[nNodes][];
		storedCorrectedPhyloTMat = new double[nSpp][nSpp];

		if (nNodes != popSizesInput.get().getValues().length) {
			throw new RuntimeException("The number of population sizes in popSizes is different from the number of nodes in the tree.");
		}
	}

	/*
	 * Does not assume tree is ultrametric
	 */
	// private int fillNLineageDistInPlace(Node aNode, Double[] allPopSizes, String[] spOrderInTMat) {
	private int fillNLineageDistInPlace(Node aNode, Double[] allPopSizes) {
		int maxNLineages = 0;
		int nodeIdx = aNode.getNr();
		Double thisNodePopSize = allPopSizes[nodeIdx];

		if (aNode.isLeaf()) {
			nLineageDistAtEnd[nodeIdx] = new double[] { 1.0 }; // always start with 1 lineage per species
			maxNLineages = 1;
			// spOrderInTMat[nodeIdx] = aNode.getID();
		}
		
		// if internal node or root
		else {
			Node leftChild = aNode.getChild(0);
			int leftIdx = leftChild.getNr();
			// int maxNLineagesLeft = fillNLineageDistInPlace(leftChild, allPopSizes, spOrderInTMat); // recursion
			int maxNLineagesLeft = fillNLineageDistInPlace(leftChild, allPopSizes); // recursion
		
			Node rightChild = aNode.getChild(1);
			int rightIdx = rightChild.getNr();
			// int maxNLineagesRight = fillNLineageDistInPlace(rightChild, allPopSizes, spOrderInTMat); // recursion
			int maxNLineagesRight = fillNLineageDistInPlace(rightChild, allPopSizes);

			// folding distributions!
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
						
						nLineageDistAtEnd[nodeIdx][kLineages-1] += probOfNLineages * 
								CoalUtils.getHeledGij(nLineages, kLineages, aNode.getLength(), thisNodePopSize);
					}
				}
			}
			
		}
		
		// System.out.println(Arrays.toString(nLineageDistAtEnd[nodeIdx]));
		return maxNLineages; 
	}
	
	/*
	 * Note that this is the expected genealogy height where the "tips" start at the root population
	 * so this is NOT the TOTAL genealogy height (i.e., this is NOT the total gene tree height, but only
	 * the height of the gene tree portion above the species tree root)
	 */
	// private double getExpGenealHeightAtRoot(Tree tree, Double[] allPopSizes, String[] spOrderInTMat) {
	private double getExpGenealHeightAtRoot(Tree tree, Double[] allPopSizes) {
		Node root = tree.getRoot();
		int rootIdx = root.getNr();
		
		double expGenealHeight = 0.0;
		// double maxNLineages = fillNLineageDistInPlace(root, allPopSizes, spOrderInTMat);
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
	
	private void fillPhyloTMatInPlace(String[] spOrderInTMat) {
		tree = treeInput.get();
		double treeHeight = tree.getRoot().getHeight();
		List<Node> allLeaves = tree.getRoot().getAllLeafNodes();
		Double[] popSizes = popSizesInput.get().getValues();
		
		/*
		 * First step: get the height of section of the gene tree that extends
		 * ABOVE the species tree root
		 */
		double expGenealHeightAtRoot = getExpGenealHeightAtRoot(tree, popSizes);

		int nSpp = tree.getLeafNodeCount();

		/*
		 * Second step: find pairwise calculations to obtain the
		 * expected times of coalescence for every pair of species
		 */
		List<Node> allNodes = tree.getRoot().getAllLeafNodes();
		for (int i=0; i<nSpp; ++i) {
			spOrderInTMat[i] = allLeaves.get(i).getID(); // this is where we populate spOrderInTMat, as opposed to in MVNUtils when coalescent correction is disabled

			for (int j=0; j<nSpp; ++j) {
				/*
				 * Diagonal entries of the matrix:
				 * Here we need to first get the height between the tip and the species tree root (treeHeight - allLeaves.get(i).getHeight(),
				 * and then add the height of the gene tree protruding ABOVE the species tree root (expGenealHeightAtRoot)
				 */
				if (i == j) {
					correctedPhyloTMat[i][j] = expGenealHeightAtRoot + (treeHeight - allLeaves.get(i).getHeight()); // all variances (diag elements) don't have to be the same
				}

				/*
				 * Off-diagonal entries of the matrix:
				 * For each pair we add the species tree height (root height) to the height of the gene tree protruding ABOVE the species tree root (expGenealHeightAtRoot)
				 */
				else {
					// we subtract the expected coalescent times from the total genealogical height (= root height + expGenealHeightAtRoot)
					correctedPhyloTMat[i][j] = tree.getRoot().getHeight() + expGenealHeightAtRoot - getExpCoalTimePair(allNodes.get(i), allNodes.get(j), popSizes);
					correctedPhyloTMat[j][i] = correctedPhyloTMat[i][j]; // assume matrix is symmetric
				}
			}
		}
	}
	
	/*
	 * Getters
	 */
	public double[][] getCorrectedPhyloTMat(String[] spOrderInTMat) {
		fillPhyloTMatInPlace(spOrderInTMat);
		
		return correctedPhyloTMat;
	}

	// caching
	@Override
	public boolean requiresRecalculation() {
		boolean dirty = false;

		if (treeInput.isDirty() || popSizesInput.isDirty()) {
			dirty = true;
		}

		return dirty;
	}

	@Override
	public void store() {
		for (int i=0; i < correctedPhyloTMat.length; ++i) {
			System.arraycopy(correctedPhyloTMat[i], 0 , storedCorrectedPhyloTMat[i], 0, correctedPhyloTMat[i].length);

			double[] nDist = new double[nLineageDistAtEnd[i].length];
			storedNLineageDistAtEnd[i] = nDist;
			for (int j=0; j < nLineageDistAtEnd[i].length; ++j) {
				System.arraycopy(nLineageDistAtEnd[i], 0, storedNLineageDistAtEnd[i], 0, nLineageDistAtEnd[i].length);
			}
		}
	}

	@Override
	public void restore() {
		double[][] array2DTmp;

		array2DTmp = correctedPhyloTMat;
		correctedPhyloTMat = storedCorrectedPhyloTMat;
		storedCorrectedPhyloTMat = array2DTmp;

		array2DTmp = nLineageDistAtEnd;
		nLineageDistAtEnd = storedNLineageDistAtEnd;
		storedNLineageDistAtEnd = array2DTmp;
	}
}
