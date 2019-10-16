package contraband;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("coalescent correction for continuous trait on tree")
public class CoalCorrection extends CalculationNode {

	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<RealParameter> popSizesInput = new Input<>("popSizes", "Temporary input containing the pop size of each node in the tree.", Validate.REQUIRED);
	
	private Tree tree;
	
	private double[][] nLineageDistAtEnd; // note that end here means earlier in "regular" (e.g., forward) time
	private double[][] correctedPhyloTMat;
	
	// array with zeros, for setting parts of nLineageDistAtEnd to 0 
	double [] nullarray;
	
	// population sizes in Double and double forms
	Double[] popSizes;
	double [] popSizesd;

	// cache for Gij calculation
    double [][][] ratesCache;
    double [][][] cisCache;

    // used for MRCA table computation
	int [][] commonAncestor; // node number of internal node that is common ancestors of pair of leafs 
	List<Integer> [] descendantsList; // memory used to prevent new ArrayList() calls

	// cache for exp(x)
	static double [] expy, //  cache for exp(0)...exp(1) in steps of 1/N
		intexp; // cache for exp(-750), exp(-749),...,exp(0)
	static int N = 1024*1024; // size of cache for exp(0)...exp(1)

	// flag to keep track of nLineageDistAtEnd cache
	// if true, we repopulate nLineageDistAtEnd entirely
    boolean hasDirt = true;

    // See with RB later
    int [] leftX, rightX, parentIdxs;
    double [] allBrLengthsInCU;

    // stored stuff (none of these seem to solve the issue)
	double[] storedPopSizesD;
	double[][] storedCorrectedPhyloMat;
	int[][] storedCommonAncestor;
	int[] storedLeftX;
	int[] storedRightX;
	int[] storedParentIdxs;
	double[] storedAllBrLengthsInCU;

	@Override
	public void initAndValidate() {
		/*
		 * exp(x) is done multiple times and takes most
		 * of the computation time; so we precompute and cache
		 * many exp(x) results
		 */
		initExp();
		
		tree = treeInput.get();
		int nSpp = tree.getLeafNodeCount();
		int nNodes = tree.getNodeCount();
		correctedPhyloTMat = new double[nSpp][nSpp];
		
		// nodesTraversed = new boolean[nNodes];
		popSizes = new Double[popSizesInput.get().getDimension()];
		popSizesd = new double[popSizes.length];
		
		nLineageDistAtEnd = new double[nNodes][];
		for (int i = 0; i < nSpp; i++) {
			nLineageDistAtEnd[i] = new double[] { 1.0 };
		}
		for (int i = nSpp; i < nNodes; i++) {
			nLineageDistAtEnd[i] = new double[nNodes];
		}
		nullarray = new double[nNodes];

		ratesCache = new double [nNodes][nNodes][];
		cisCache = new double [nNodes][nNodes][];
		
		commonAncestor = new int[nSpp][nSpp];
		descendantsList = new List[nNodes];
		// every node will have an entry in descendantsList
		for (int i = 0; i < nNodes; i++) {
			descendantsList[i] = new ArrayList<>();
		}
		for (int i = 0; i < nSpp; i++) {
			descendantsList[i].add(i); // we pre-initialize the tip nodes with their own node numbers
		}		
//		storedNLineageDistAtEnd = new double[nNodes][];
//		storedCorrectedPhyloTMat = new double[nSpp][nSpp];

		leftX = new int[nNodes];
		rightX = new int[nNodes];
		allBrLengthsInCU = new double[nNodes]; // CU = coalescent units
		parentIdxs = new int[nNodes];
		
		if (nNodes != popSizesInput.get().getValues().length) {
			throw new RuntimeException("The number of population sizes in popSizes is different from the number of nodes in the tree.");
		}

		// stored stuff
		storedPopSizesD = new double[popSizes.length];
		storedCorrectedPhyloMat = new double[nSpp][nSpp];
		storedCommonAncestor = new int[nSpp][nSpp];
		storedLeftX = new int[nNodes];
		storedRightX = new int[nNodes];
		storedParentIdxs = new int[nNodes];
		storedAllBrLengthsInCU = new double[nNodes];
	}

	/*
	 * Does not assume tree is ultrametric
	 */
	// private int fillNLineageDistInPlace(Node aNode, Double[] allPopSizes, String[] spOrderInTMat) {
	private int fillNLineageDistInPlace(Node aNode) {
		if (aNode.isLeaf()) {
			// if node is sampled ancestor, it contributes no lineages
			if (aNode.isDirectAncestor()) {
				return 0;
			}

			/*
			 * RB code: This allows caching of nLineageDistAtEnd
			 *
			 * Here, we return the number of lineages a leaf contributes (1),
			 * but also information on the dirtiness of the leaf (the negative sign)
			 *
			 * Importantly, we always need to flip the sign of the number of lineages
			 * coming into an internal node from a leaf
			 */
			int parentIdx = aNode.getParent().getNr();
			// if my parent changed (if I was moved to a new place)
			if (parentIdxs[aNode.getNr()] != parentIdx) {
				parentIdxs[aNode.getNr()] = parentIdx;
				return -1; // ask RB how this works for SA's, as those contribute no lineages
			}
            return 1;

			// original code
			// nLineageDistAtEnd[nodeIdx] = new double[] { 1.0 }; // always start with 1 lineage per species
			// maxNLineages = 1;
			// spOrderInTMat[nodeIdx] = aNode.getID();
		}
		
		// if internal node or root
		else {
            boolean isDirty = false;
            
            int maxNLineages = 0;
			int nodeIdx = aNode.getNr();
			double thisNodePopSize = popSizesd[nodeIdx];

			Node leftChild = aNode.getChild(0);
			int leftIdx = leftChild.getNr();

			// original code
			// int maxNLineagesLeft = fillNLineageDistInPlace(leftChild, allPopSizes, spOrderInTMat); // recursion

			// RB code
			int maxNLineagesLeft = fillNLineageDistInPlace(leftChild); // recursion
			if (maxNLineagesLeft < 0) {
				/*
				 * If sign was negative, it means left child has been operated on,
				 * or is in the path between a node that has been operated on and the root
				 */
				isDirty = true;
                maxNLineagesLeft = -maxNLineagesLeft; // flipping sign back
            }

			Node rightChild = aNode.getChild(1);
			int rightIdx = rightChild.getNr();

			// original code
			// int maxNLineagesRight = fillNLineageDistInPlace(rightChild, allPopSizes, spOrderInTMat); // recursion

			int maxNLineagesRight = fillNLineageDistInPlace(rightChild);
            if (maxNLineagesRight < 0) {
                // same as above
            	isDirty = true;
                maxNLineagesRight = -maxNLineagesRight;
            }

			int parent = aNode.isRoot() ? -1 : aNode.getParent().getNr();
            double brLengthCU = -aNode.getLength() / thisNodePopSize;
            if (parentIdxs[nodeIdx] != parent ||
            	leftX[nodeIdx] != maxNLineagesLeft || 
            	rightX[nodeIdx] != maxNLineagesRight || 
            	allBrLengthsInCU[nodeIdx] != brLengthCU) {
            	isDirty = true;
            }
            parentIdxs[nodeIdx] = parent;
            leftX[nodeIdx] = maxNLineagesLeft; 
            rightX[nodeIdx] = maxNLineagesRight; 
            allBrLengthsInCU[nodeIdx] = brLengthCU;

			// folding distributions below
			// original code
			// nLineageDistAtEnd[nodeIdx] = new double[maxNLineages];
			maxNLineages = maxNLineagesLeft + maxNLineagesRight;
			if (true || isDirty) { // ignoring the caching of nLineageDistAtEnd while fixing sampled ancestors part
				double [] left = nLineageDistAtEnd[leftIdx];
				double [] right = nLineageDistAtEnd[rightIdx];
				double [] nodel = nLineageDistAtEnd[nodeIdx];
			
				// fill nodel with zeros
				System.arraycopy(nullarray, 0, nodel, 0, maxNLineages);
				// Arrays.fill(nodel, 0, nodel.length-1, 0.0); // ask Remco later if this is good

				// see comments for when neither children was sampled ancestors below
				// for an explanation of the code

				// if left child was sampled ancestor
				if (maxNLineagesLeft == 0) {
					for (int nLineages=1; nLineages <= maxNLineages; ++nLineages) {
						double probOfNLineages = 0.0;

						for (int nRight=1; nRight <= maxNLineagesRight; ++nRight) {
							probOfNLineages += right[nRight-1];
						}

						if (aNode.isRoot()) {
							nodel[nLineages-1] = probOfNLineages;
						}

						else {
							for (int kLineages=1; kLineages <= nLineages; ++kLineages) {
								nodel[kLineages-1] += probOfNLineages *
										getHeledGij(nLineages, kLineages, brLengthCU);
							}
						}
					}
				}

				// if right child was sampled ancestor
				else if (maxNLineagesRight == 0) {
					for (int nLineages=1; nLineages <= maxNLineages; ++nLineages) {
						double probOfNLineages = 0.0;

						for (int nLeft=1; nLeft <= maxNLineagesLeft; ++nLeft) {
							probOfNLineages += left[nLeft - 1];
						}

						if (aNode.isRoot()) {
							nodel[nLineages - 1] = probOfNLineages;
						}

						else {
							for (int kLineages = 1; kLineages <= nLineages; ++kLineages) {
								nodel[kLineages - 1] += probOfNLineages *
										getHeledGij(nLineages, kLineages, brLengthCU);
							}
						}
					}
				}

				// neither children was sampled ancestor
				else {
					/*
					 * nLineages starts at 2 because if neither children are sampled ancestors,
					 * we have a minimum of 1 lineage coming from either side (=2)
					 */
					for (int nLineages=2; nLineages <= maxNLineages; ++nLineages) {
						/*
						 * We are going backward in time
						 *
						 * We start at nLineages with probOfNLineages,
						 * which we initialize at 0.0 below, but update by considering all combinations
						 * of left and right lineages that summed up to nLineages
						 */
						double probOfNLineages = 0.0;
						for (int nLeft=1; nLeft <= maxNLineagesLeft; ++nLeft) {
							int nRight = nLineages - nLeft; // this allows us to consider different combinations, as we're always setting nRight to be the reciprocal of nLeft
						
							if (nRight > 0 && nRight <= maxNLineagesRight) {
								/*
								 * Once we determine a combination of lineages from left and right,
								 * we now get the probability of that combination, and add it to (the
								 * running sum) probOfNLineages
								 */
								probOfNLineages += left[nLeft-1] * right[nRight-1]; // -1 is offset
							}
						}
					
						// at the root, we stop at just the probOfNLineages at the start, no Tavare coeffient
						if (aNode.isRoot()) {
							nodel[nLineages-1] = probOfNLineages;
						}

						/*
						 * At other internal nodes, we have nLineages with probOfNLineages at the start, and then
						 * these nLineages might coalesce to kLineages according to Tavare's coefficients
						 */
						else {
							for (int kLineages=1; kLineages <= nLineages; ++kLineages) {

								nodel[kLineages-1] += probOfNLineages *
										getHeledGij(nLineages, kLineages, brLengthCU);
							}
						}
					}
			 	}
			}

			/*
			 * RB code: as we do for leaf nodes, we return maxNLineages for this internal node,
			 * but add a negative sign if this internal node was operated on, or if its on a path from
			 * a node operated on and the root
			 *
			 * As before, this sign has to be flipped when we recur back to this internal node's parent
			 */
			if (isDirty) {
				return -maxNLineages;
			} else {
				return maxNLineages;
			}
		}
		
		// System.out.println(Arrays.toString(nLineageDistAtEnd[nodeIdx]));
	}
	
	/*
	 * Note that this is the expected genealogy height where the "tips" start at the root population
	 * so this is NOT the TOTAL genealogy height (i.e., this is NOT the total gene tree height, but only
	 * the height of the gene tree portion above the species tree root)
	 */
	// private double getExpGenealHeightAtRoot(Tree tree, Double[] allPopSizes, String[] spOrderInTMat) {
	private double getExpGenealHeightAtRoot(Tree tree) {
		Node root = tree.getRoot();
		int rootIdx = root.getNr();
		
		double expGenealHeight = 0.0;
		// double maxNLineages = fillNLineageDistInPlace(root, allPopSizes, spOrderInTMat);
		double maxNLineages = Math.abs(fillNLineageDistInPlace(root)); // call Math.abs because it could have a negative sign coming from the caching of nLineageDistAtEnd
		hasDirt = false;
		
		for (int kLineages=1; kLineages <= maxNLineages; ++kLineages) {
			double probOfKLineages = nLineageDistAtEnd[rootIdx][kLineages-1];
			expGenealHeight += probOfKLineages * CoalUtils.getMeanRootHeight(kLineages, popSizesd[rootIdx]);
		}
		
		return expGenealHeight;
	}

	private double getExpCoalTimePair(Node mrcaInSpTree, Double[] allPopSizes) {
		int k = mrcaInSpTree.getNr();
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
			} else {
				e = (1 - (branchLength * invPopSize + 1) * Math.exp(-branchLength * invPopSize)) / invPopSize;
				probNoCoal = getHeledGij(2, 2, -branchLength/popSize);
			}
			
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
		//List<Node> allLeaves = tree.getRoot().getAllLeafNodes();
		
		 popSizesInput.get().getValues(popSizes);
		 for (int i = 0; i < popSizes.length; i++) {
			 popSizesd[i] = popSizes[i];
		 }
		 
		
		/*
		 * First step: get the height of section of the gene tree that extends
		 * ABOVE the species tree root
		 */
		double expGenealHeightAtRoot = getExpGenealHeightAtRoot(tree);

		int nSpp = tree.getLeafNodeCount();
		calcCommonAncestors(tree.getRoot());

		/*
		 * Second step: find pairwise calculations to obtain the
		 * expected times of coalescence for every pair of species
		 */
		//List<Node> allNodes = tree.getRoot().getAllLeafNodes();
		Node [] nodes = tree.getNodesAsArray();
		double [] expCoalTimePair = new double[nodes.length];
		for (int i = tree.getLeafNodeCount(); i < nodes.length; i++) {
			expCoalTimePair[i] = getExpCoalTimePair(nodes[i], popSizes);
		}
		
		
		for (int i=0; i<nSpp; ++i) {
			spOrderInTMat[i] = nodes[i].getID(); // this is where we populate spOrderInTMat, as opposed to in MVNUtils when coalescent correction is disabled

			for (int j=i; j<nSpp; ++j) {
				/*
				 * Diagonal entries of the matrix:
				 * Here we need to first get the height between the tip and the species tree root (treeHeight - allLeaves.get(i).getHeight(),
				 * and then add the height of the gene tree protruding ABOVE the species tree root (expGenealHeightAtRoot)
				 */
				if (i == j) {
					correctedPhyloTMat[i][j] = expGenealHeightAtRoot + (treeHeight - nodes[i].getHeight()); // all variances (diag elements) don't have to be the same
				}

				/*
				 * Off-diagonal entries of the matrix:
				 * For each pair we add the species tree height (root height) to the height of the gene tree protruding ABOVE the species tree root (expGenealHeightAtRoot)
				 */
				else {
					// we subtract the expected coalescent times from the total genealogical height (= root height + expGenealHeightAtRoot)
					// correctedPhyloTMat[i][j] = tree.getRoot().getHeight() + expGenealHeightAtRoot - getExpCoalTimePair(nodes[i], nodes[j], popSizes);
					correctedPhyloTMat[i][j] = treeHeight + expGenealHeightAtRoot - expCoalTimePair[commonAncestor[i][j]];
					correctedPhyloTMat[j][i] = correctedPhyloTMat[i][j]; // assume matrix is symmetric
				}
			}
		}
	}
	
	/*
	 * Getters
	 */
	public double[][] getCorrectedPhyloTMat(String[] spOrderInTMat) {
		
		if (tree.getRoot().getNr() != tree.getNodeCount()-1) {
			int h = 3;
			h++;
		}
		fillPhyloTMatInPlace(spOrderInTMat);
		
		return correctedPhyloTMat;
	}


	
	
	class X {
		Integer from, to; Double scale;
		X(int from, int to, double scale) {
			this.from = from; this.to = to; this.scale = scale;
		}
	}
	
	Map<X, Double> GijCache = new HashMap<>();
    
	private double getHeledGij(int nFrom, int nTo, double scale) {
//		double scale = -t/pop;
//		X key = new X(nFrom, nTo, scale);
//		if (GijCache.containsKey(key)) {
//			return GijCache.get(key);
//		}
		
		
		
		if (ratesCache[nFrom][nTo] == null) {
		
			// (n-k) rates, backwards
			double[] rates = new double[nFrom - nTo + 1];
			int counter = 0;
			for (int i=nTo; i < (nFrom+1); ++i) {
				rates[counter] = (i*(i-1.0))/2.0;
				counter++;
			}
			
			double[] cis = new double[nFrom - nTo + 1];
			cis[0] = 1.0;
			int cisCounter = 1;
			for (int i=1; i < (nFrom - nTo + 1); ++i) {
				double rate = rates[i];
				
				for (int j=0; j < cisCounter; ++j) {
					cis[j] *= (rate / (rate - rates[j]));
				}
	
				// cis[cisCounter] = -Arrays.stream(cis).sum();
				double sum = 0;
				for (double d : cis) {
					sum += d;
				}
				cis[cisCounter] = -sum;
				
				cisCounter++;
			}
			ratesCache[nFrom][nTo] = rates;
			cisCache[nFrom][nTo] = cis;
		}
		
		double [] rates = ratesCache[nFrom][nTo];
		double [] cis = cisCache[nFrom][nTo];
		
		double prob = 0.0;
		for (int i=0; i < rates.length; ++i) {
			//prob += cis[i] * FastMath. exp(rates[i] * scale);
			prob += cis[i] * exp(rates[i] * scale);
		}
		
//		System.out.println("Printing from n to k lineages: n=" + nFrom + " k=" + nTo + " d=" + pop + " over all interval, t=" + t);
//		System.out.println("n to k lineages result=" + prob);

//		GijCache.put(key, prob);
		return prob;
	}

	/*
	 * RB code: When populating phyloTMat pairwise, we need to grab
	 * the common ancestor of two different nodes several times, which requires
	 * multiple tree traversals if this is done by a function call
	 *
	 * Here, we precompute all the common ancestors for all pairs of species and
	 * store them in the 2D-int array "commonAncestor". So there is a single tree
	 * traversal done. We then just query common ancestors from "commonAncestor"
	 */
    protected List<Integer> calcCommonAncestors(Node node) {
		int nodeIdx = node.getNr();

		if (node.isLeaf()) {
			// had already been populated in initAndValidate
			// returns leaf index
    		return descendantsList[nodeIdx];
    	}

		// if internal node or root
    	else {
    		// first we list all children
    		List<Integer> left = calcCommonAncestors(node.getLeft()); // recursive step
    		List<Integer> right = calcCommonAncestors(node.getRight()); // recursive step

			// now we iterate over all pairs of left and right tip children
    		for (int j : left) {
    			for (int k : right) {
    				// populating the actual MRCA table
    				commonAncestor[j][k] = nodeIdx; // the current node is the actual MRCA of j and k
    				commonAncestor[k][j] = nodeIdx; // the current node is the actual MRCA of k and j
    			}
    		}

    		List<Integer> myDescList = descendantsList[nodeIdx]; // in initAndValidate we initialized it to an empty ArrayList
    		myDescList.clear();
    		myDescList.addAll(left); // adding left descendants
    		myDescList.addAll(right); // adding right descendants
    		return myDescList;
    	}
    }
	
    // exponential calculation
	private void initExp() {
		expy = new double[N];
		for (int i = 0; i < N; i++) {
			expy[i] = Math.exp(-i/(N - 1.0));
		}
		intexp = new double[750];
		for (int i = 0; i < 750; i++) {
			intexp[i] = Math.exp(-i);
		}
	}
	
    // exponential calculation
	// assumes val < 0
	private double exp(double val) {
		if (val < -746) {
			return 0;
		}
		int v = (int) val;
		double rest = - val % 1.0;
		int i = (int)(rest * N);
		return expy[i] * intexp[-v];
	}


	// caching
	@Override
	public boolean requiresRecalculation() {
		// hasDirt = false;
		// hasDirt = true;

		if (popSizesInput.isDirty()) {
			hasDirt = true;
			return true;
		}
	    if (treeInput.get().somethingIsDirty()) {
	    	return true;
	    }
		
		return false;
	}
	
	@Override
	protected void store() {
		// hasDirt = true;
		int nSpp = tree.getLeafNodeCount();
		for (int i=0; i<nSpp; ++i) {
			for (int j=0; j<nSpp; ++j) {
				storedCorrectedPhyloMat[i][j] = correctedPhyloTMat[i][j];
				storedCommonAncestor[i][j] = commonAncestor[i][j];
			}
		}

		for (int i=0; i<tree.getNodeCount(); ++i) {
			storedLeftX[i] = leftX[i];
			storedRightX[i] = rightX[i];
			storedParentIdxs[i] = parentIdxs[i];
			storedAllBrLengthsInCU[i] = allBrLengthsInCU[i];
		}

		for (int i=0; i<popSizesInput.get().getDimension(); ++i) {
			storedPopSizesD[i] = popSizesd[i];
		}

		super.store();
	}
	
	@Override
	protected void restore() {
    	double[] doubleArrTmp;
		double[][] double2DArrTmp;
		int[] intArrTmp;
		int[][] int2DArrTmp;

    	doubleArrTmp = popSizesd;
    	popSizesd = storedPopSizesD;
    	storedPopSizesD = doubleArrTmp;

    	double2DArrTmp = correctedPhyloTMat;
    	correctedPhyloTMat = storedCorrectedPhyloMat;
    	storedCorrectedPhyloMat = double2DArrTmp;

    	intArrTmp = leftX;
    	leftX = storedLeftX;
    	storedLeftX = intArrTmp;

		intArrTmp = rightX;
		rightX = storedRightX;
		storedRightX = intArrTmp;

		intArrTmp = parentIdxs;
		parentIdxs = storedParentIdxs;
		storedParentIdxs = intArrTmp;

		doubleArrTmp = allBrLengthsInCU;
		allBrLengthsInCU = storedAllBrLengthsInCU;
		storedAllBrLengthsInCU = doubleArrTmp;

		int2DArrTmp = commonAncestor;
		commonAncestor = storedCommonAncestor;
		storedCommonAncestor = int2DArrTmp;

		hasDirt = true;
		super.restore();
	}
}
