package contraband;


import java.util.ArrayList;
import java.util.Arrays;
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

    
	int [][] commonAncestor; // node number of internal node that is common ancestors of pair of leafs 
	List<Integer> [] commonAncestorList; // memory used to prevent new ArrayList() calls

	// cache for exp(x)
	static double [] expy, //  cache for exp(0)...exp(1) in steps of 1/N
		intexp; // cache for exp(-750), exp(-749),...,exp(0)
	static int N = 1024*1024; // size of cache for exp(0)...exp(1)


	// flag to keep track of cache
    boolean hasDirt = true;

	@Override
	public void initAndValidate() {
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
		commonAncestorList = new List[nNodes];
		for (int i = 0; i < nNodes; i++) {
			commonAncestorList[i] = new ArrayList<>();
		}
		for (int i = 0; i < nSpp; i++) {
			commonAncestorList[i].add(i);
		}		
//		storedNLineageDistAtEnd = new double[nNodes][];
//		storedCorrectedPhyloTMat = new double[nSpp][nSpp];

		if (nNodes != popSizesInput.get().getValues().length) {
			throw new RuntimeException("The number of population sizes in popSizes is different from the number of nodes in the tree.");
		}
	}

	/*
	 * Does not assume tree is ultrametric
	 */
	// private int fillNLineageDistInPlace(Node aNode, Double[] allPopSizes, String[] spOrderInTMat) {
	private int fillNLineageDistInPlace(Node aNode) {

		if (aNode.isLeaf()) {
			if (aNode.isDirectAncestor()) {
				return 0;
			}
			// nLineageDistAtEnd[nodeIdx] = new double[] { 1.0 }; // always start with 1 lineage per species
			// maxNLineages = 1;
			// spOrderInTMat[nodeIdx] = aNode.getID();
			if (hasDirt || (aNode.isDirty() != Tree.IS_CLEAN)) {
                return -1;
			} else {
                return 1;
			}		
		}
		
		// if internal node or root
		else {
            boolean isDirty = hasDirt || (aNode.isDirty() != Tree.IS_CLEAN);
            
            int maxNLineages = 0;
			int nodeIdx = aNode.getNr();
			double thisNodePopSize = popSizesd[nodeIdx];

			Node leftChild = aNode.getChild(0);
			int leftIdx = leftChild.getNr();
			// int maxNLineagesLeft = fillNLineageDistInPlace(leftChild, allPopSizes, spOrderInTMat); // recursion
			int maxNLineagesLeft = fillNLineageDistInPlace(leftChild); // recursion
            if (maxNLineagesLeft < 0) {
                isDirty = true;
                maxNLineagesLeft = -maxNLineagesLeft;
            }

			Node rightChild = aNode.getChild(1);
			int rightIdx = rightChild.getNr();
			// int maxNLineagesRight = fillNLineageDistInPlace(rightChild, allPopSizes, spOrderInTMat); // recursion
			int maxNLineagesRight = fillNLineageDistInPlace(rightChild);
            if (maxNLineagesRight < 0) {
                isDirty = true;
                maxNLineagesRight = -maxNLineagesRight;
            }

			// folding distributions!
			maxNLineages = maxNLineagesLeft + maxNLineagesRight;
			// nLineageDistAtEnd[nodeIdx] = new double[maxNLineages];
			
			if (true || isDirty) {
			double [] left = nLineageDistAtEnd[leftIdx];
			double [] right = nLineageDistAtEnd[rightIdx];
			double [] nodel = nLineageDistAtEnd[nodeIdx];
			
			// fill nodel with zeros
			System.arraycopy(nullarray, 0, nodel, 0, maxNLineages);
			//Arrays.fill(nodel, 0, 0, maxNLineages);
			
			if (maxNLineagesLeft == 0) {
				
			} else if (maxNLineagesRight == 0) {
			
			} else {
				for (int nLineages=2; nLineages <= maxNLineages; ++nLineages) {
					double probOfNLineages = 0.0; // going backward in time, we start at nLineages with this prob
					
					for (int nLeft=1; nLeft <= maxNLineagesLeft; ++nLeft) {
						int nRight = nLineages - nLeft;
						
						if (nRight > 0 && nRight <= maxNLineagesRight) {
							// we're iterating over all combinations of nLeft and nRight that comprise values from 2 to maxNLineages
							probOfNLineages += left[nLeft-1] * right[nRight-1]; // -1 is offset
						}
					}
					
					// at the root, we stop at just the probOfNLineages at the start, no Tavare coeffient
					if (aNode.isRoot()) {
						nodel[nLineages-1] = probOfNLineages;
					}
					// at other internal nodes, we have nLineages with probOfNLineages at the start, and then
					// these nLineages might coalesce to kLineages according to Tavare's coefficients
					else {
						for (int kLineages=1; kLineages <= nLineages; ++kLineages) {
							
							nodel[kLineages-1] += probOfNLineages * 
									getHeledGij(nLineages, kLineages, aNode.getLength(), thisNodePopSize);
						}
					}
				}
			}
			}
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
		double maxNLineages = Math.abs(fillNLineageDistInPlace(root));
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
				probNoCoal = getHeledGij(2, 2, branchLength, popSize);
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
    
	private double getHeledGij(int nFrom, int nTo, double t, double pop) {
		double scale = -t/pop;
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

    protected List<Integer> calcCommonAncestors(Node node) {
		int i = node.getNr();
    	if (node.isLeaf()) {
    		return commonAncestorList[i];
    	} else {
    		List<Integer> left = calcCommonAncestors(node.getLeft());
    		List<Integer> right = calcCommonAncestors(node.getRight());
    		for (int j : left) {
    			for (int k : right) {
    				commonAncestor[j][k] = i;
    				commonAncestor[k][j] = i;
    			}
    		}
    		List<Integer> list = commonAncestorList[i];
    		list.clear();
    		list.addAll(left);
    		list.addAll(right);
    		return list;
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
		hasDirt = false;
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
		super.store();
	}
	
	@Override
	protected void restore() {
		hasDirt = true;
		super.restore();
	}
}
