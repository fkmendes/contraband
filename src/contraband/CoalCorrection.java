package contraband;

import java.util.List;

import org.apache.commons.math3.util.FastMath;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class CoalCorrection extends CalculationNode {

	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<RealParameter> popSizesInput = new Input<>("popSizes", "Temporary input containing the pop size of each node in the tree.", Validate.REQUIRED);
	
	private Tree tree;
	
	private double[][] nLineageDistAtEnd; // note that end here means earlier in "regular" (e.g., forward) time
	private double[][] correctedPhyloTMat;
	
	double [] nullarray;

	// stored stuff
//	private double[][] storedCorrectedPhyloTMat;
//	private double[][] storedNLineageDistAtEnd;

	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		int nSpp = tree.getLeafNodeCount();
		int nNodes = tree.getNodeCount();
		correctedPhyloTMat = new double[nSpp][nSpp];
		
		nodesTraversed = new boolean[nNodes];
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
			// nLineageDistAtEnd[nodeIdx] = new double[] { 1.0 }; // always start with 1 lineage per species
			// maxNLineages = 1;
			// spOrderInTMat[nodeIdx] = aNode.getID();
			return 1;
		}
		
		// if internal node or root
		else {
			int maxNLineages = 0;
			int nodeIdx = aNode.getNr();
			double thisNodePopSize = popSizesd[nodeIdx];

			Node leftChild = aNode.getChild(0);
			int leftIdx = leftChild.getNr();
			// int maxNLineagesLeft = fillNLineageDistInPlace(leftChild, allPopSizes, spOrderInTMat); // recursion
			int maxNLineagesLeft = fillNLineageDistInPlace(leftChild); // recursion
		
			Node rightChild = aNode.getChild(1);
			int rightIdx = rightChild.getNr();
			// int maxNLineagesRight = fillNLineageDistInPlace(rightChild, allPopSizes, spOrderInTMat); // recursion
			int maxNLineagesRight = fillNLineageDistInPlace(rightChild);

			// folding distributions!
			maxNLineages = maxNLineagesLeft + maxNLineagesRight;
			// nLineageDistAtEnd[nodeIdx] = new double[maxNLineages];
			
			double [] left = nLineageDistAtEnd[leftIdx];
			double [] right = nLineageDistAtEnd[rightIdx];
			double [] nodel = nLineageDistAtEnd[nodeIdx];
			
			// fill nodel with zeros
			System.arraycopy(nullarray, 0, nodel, 0, maxNLineages);
			
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
			return maxNLineages; 
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
		double maxNLineages = fillNLineageDistInPlace(root);
		
		for (int kLineages=1; kLineages <= maxNLineages; ++kLineages) {
			double probOfKLineages = nLineageDistAtEnd[rootIdx][kLineages-1];
			expGenealHeight += probOfKLineages * CoalUtils.getMeanRootHeight(kLineages, popSizesd[rootIdx]);
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
//		Set<String> pairOfSpNames = new HashSet<String>();
//		pairOfSpNames.add(node1.getID());
//		pairOfSpNames.add(node2.getID());
//		Node mrcaInSpTree = TreeUtils.getCommonAncestorNode(tree, pairOfSpNames); // mrca of two nodes in species tree
		
		java.util.Arrays.fill(nodesTraversed, false);
		Node mrcaInSpTree = getCommonAncestor(node1, node2);
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
	
	Double[] popSizes;
	double [] popSizesd;

	private void fillPhyloTMatInPlace(String[] spOrderInTMat) {
		tree = treeInput.get();
		double treeHeight = tree.getRoot().getHeight();
		List<Node> allLeaves = tree.getRoot().getAllLeafNodes();
		
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

//	// caching
//	@Override
//	public boolean requiresRecalculation() {
//		boolean dirty = false;
//
//		if (treeInput.isDirty() || popSizesInput.isDirty()) {
//			dirty = true;
//		}
//
//		return dirty;
//	}
//
//	@Override
//	public void store() {
//		for (int i=0; i < correctedPhyloTMat.length; ++i) {
//			System.arraycopy(correctedPhyloTMat[i], 0 , storedCorrectedPhyloTMat[i], 0, correctedPhyloTMat[i].length);
//
//			double[] nDist = new double[nLineageDistAtEnd[i].length];
//			storedNLineageDistAtEnd[i] = nDist;
//			for (int j=0; j < nLineageDistAtEnd[i].length; ++j) {
//				System.arraycopy(nLineageDistAtEnd[i], 0, storedNLineageDistAtEnd[i], 0, nLineageDistAtEnd[i].length);
//			}
//		}
//	}
//
//	@Override
//	public void restore() {
//		double[][] array2DTmp;
//
//		array2DTmp = correctedPhyloTMat;
//		correctedPhyloTMat = storedCorrectedPhyloTMat;
//		storedCorrectedPhyloTMat = array2DTmp;
//
//		array2DTmp = nLineageDistAtEnd;
//		nLineageDistAtEnd = storedNLineageDistAtEnd;
//		storedNLineageDistAtEnd = array2DTmp;
//	}
	

	// used by getCommonAncestor
	boolean [] nodesTraversed;
	
    protected Node getCommonAncestor(Node n1, Node n2) {
        // assert n1.getTree() == n2.getTree();
        if( ! nodesTraversed[n1.getNr()] ) {
            nodesTraversed[n1.getNr()] = true;
        }
        if( ! nodesTraversed[n2.getNr()] ) {
            nodesTraversed[n2.getNr()] = true;
        }
        while (n1 != n2) {
	        double h1 = n1.getHeight();
	        double h2 = n2.getHeight();
	        if ( h1 < h2 ) {
	            n1 = n1.getParent();
	            if( ! nodesTraversed[n1.getNr()] ) {
	                nodesTraversed[n1.getNr()] = true;
	            }
	        } else if( h2 < h1 ) {
	            n2 = n2.getParent();
	            if( ! nodesTraversed[n2.getNr()] ) {
	                nodesTraversed[n2.getNr()] = true;
	            }
	        } else {
	            //zero length branches hell
	            Node n;
	            double b1 = n1.getLength();
	            double b2 = n2.getLength();
	            if( b1 > 0 ) {
	                n = n2;
	            } else { // b1 == 0
	                if( b2 > 0 ) {
	                    n = n1;
	                } else {
	                    // both 0
	                    n = n1;
	                    while( n != null && n != n2 ) {
	                        n = n.getParent();
	                    }
	                    if( n == n2 ) {
	                        // n2 is an ancestor of n1
	                        n = n1;
	                    } else {
	                        // always safe to advance n2
	                        n = n2;
	                    }
	                }
	            }
	            if( n == n1 ) {
                    n = n1 = n.getParent();
                } else {
                    n = n2 = n.getParent();
                }
	            if( ! nodesTraversed[n.getNr()] ) {
	                nodesTraversed[n.getNr()] = true;
	            } 
	        }
        }
        return n1;
    }


    double [][][] ratesCache;
    double [][][] cisCache;
    
	private double getHeledGij(int nFrom, int nTo, double t, double pop) {
		
		if (ratesCache[nFrom][nTo] == null) {
		
			// nFrom = n
			// nTo = k
			
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
		double scale = -t/pop;
		for (int i=0; i < rates.length; ++i) {
			prob += cis[i] * FastMath. exp(rates[i] * scale);
			// prob += cis[i] * exp3(rates[i] * scale);
		}
		
//		System.out.println("Printing from n to k lineages: n=" + nFrom + " k=" + nTo + " d=" + pop + " over all interval, t=" + t);
//		System.out.println("n to k lineages result=" + prob);
				
		return prob;
	}

	
	double exp3(double x) {
	  x = 1.0 + x / 1024;
	  x *= x; x *= x; x *= x; x *= x;
	  x *= x; x *= x; x *= x; x *= x;
	  x *= x; x *= x;
	  return x;
	}
	
	double exp2(double x) {
		// for negative argument, use identity e^-x =  1/e^x
        boolean isNegative = false;
        if (x < 0) {
            isNegative = true;
            x = -x;
        }

        // compute e^x assuming x >= 0
        double term = 1.0;
        double sum = 0.0;
        for (int n = 1; sum != sum + term; n++) {
            sum += term;
            term *= x/n;
        }

        // print results
        if (isNegative)
            sum = 1.0 / sum;
        return sum;
	}
	
	
	public static double exp(double val) {
	    final long tmp = (long) (1512775 * val + (1072693248 - 60801));
	    return Double.longBitsToDouble(tmp << 32);
	}
}
