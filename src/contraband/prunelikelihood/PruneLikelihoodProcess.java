package contraband.prunelikelihood;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import java.util.List;
import java.util.Random;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.utils.PruneLikelihoodUtils;


public abstract class PruneLikelihoodProcess extends Distribution {
    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Input.Validate.REQUIRED);
    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel", "the rate or optimum on each branch");
    final public Input<NodeMath> nodeMathInput = new Input<>("nodeMath","Node information that will be used in PCM likelihood calculation.", Input.Validate.REQUIRED);
    final public Input<RealParameter> traitsValuesInput = new Input<>("traits","Trait values at tips.", Input.Validate.REQUIRED);
    //final public Input<KeyRealParameter> traitsValuesInput = new Input<>("traits","Trait values at tips.", Input.Validate.REQUIRED);


    private Tree tree;
    private int nTraits;
    private int nSpecies;
    private int nSpeciesWithData;

    private RealParameter traitsValues;

    private BranchRateModel.Base branchRateModel;

    private NodeMath nodeMath;

    private double[] traitValuesArr;

    private boolean popSE;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // get the tree
        tree = treeInput.get();

        // get clock model
        branchRateModel = branchRateModelInput.get();

        // get the trait values for tips
        // and make a list of real vectors
        traitsValues = traitsValuesInput.get();
        nTraits = traitsValues.getMinorDimension1();
        nSpeciesWithData = traitsValues.getMinorDimension2();
        nSpecies = tree.getLeafNodeCount();
        traitValuesArr = new double[nSpecies * nTraits];

        // check input
        if (nodeMathInput.get() == null) {
            throw new RuntimeException("PruneLikelihoodProcess::NodeMath is required for pmc likelihood.");
        }
        nodeMath = nodeMathInput.get();

        // populate the trait values in an array
        // each species has nTraits in the array
        //PruneLikelihoodUtils.populateTraitValuesArr(traitsValues, tree, nTraits, traitValuesArr);
        PruneLikelihoodUtils.populateTraitValuesArr(traitsValues, tree, nodeMath, nTraits, traitValuesArr);
    }

    protected void populateLogP() {

        // reject the tree with a sampled ancestor being the child of root
        if(tree.getRoot().getChild(0).isDirectAncestor() || tree.getRoot().getChild(1).isDirectAncestor()) {
            logP =  Double.NEGATIVE_INFINITY;
            return;
        }
        // prune the tree by starting from the root
        nodeMath.setLikelihoodForSampledAncestors(0.0);

        // to obtain the node information about missing data
        nodeMath.initializeNodeStatArrays();
        prune(tree.getRoot());

        // if using shrinkage method, 'traitValuesArr' is 'transformedTraitValues'
        // otherwise, it is original trait values.
        pruneNode(tree.getRoot(), nTraits, traitValuesArr, branchRateModel, nodeMath, popSE);

        // if at some internal node, (aMat + lMat) is singular or -2 * (aMat + lMat) is singular,
        // when calculating (aMat + lMat).inverse and det[-2 * (aMat + lMat)],
        // reject the state
        if (nodeMath.isSingularMatrix()) {
            logP = Double.NEGATIVE_INFINITY;
            return;
        }

        // get L, m and r at the root
        int rootIdx = tree.getRoot().getNr();
        double l0 = nodeMath.getLForNode(rootIdx);
        double[] m0 = nodeMath.getMVecForNode(rootIdx);
        double r0 = nodeMath.getRForNode(rootIdx);

        // get the root values
        nodeMath.populateRootValuesVec(rootIdx);

        // calculate likelihood in log space
        logP = calculateLikelihood(nodeMath, l0, m0, r0, rootIdx);
    }

    // getters
    public int getNTraits (){ return nTraits; }

    public int getNSpecies (){ return nSpecies; }

    public double getLogP () { return logP; }

    public NodeMath getNodeMath () { return nodeMath; }

    public int getRootIndex () { return tree.getRoot().getNr(); }

    public int getNumberOfSpeciesWithData() {return nSpeciesWithData; }

    // setters
    public void setPopSE (boolean value) { popSE = value; }

    public void setTraitValuesArr (double[] values) {
        traitValuesArr = values;
    }

    public void pruneNode(Node node, int nTraits, double[] traitValuesArr,
                           BranchRateModel.Base pcmc, NodeMath nodeMath, boolean popSE) {

        int thisNodeIdx = node.getNr();

        // initialize before sum
        double thisNodeL = 0.0;
        double [] thisNodeMVec = nodeMath.getInitMVec().clone();
        double thisNodeR = 0.0;

        List<Node> children = node.getChildren();

        for (Node child : children) {

            int childIdx = child.getNr();

            if(nodeMath.isSpeciesToIgnore(childIdx)){
                continue;
            }

            double branchLength = child.getLength() * pcmc.getRateForBranch(child);
            nodeMath.setVarianceForTip(childIdx, branchLength);

            // (1) child is a normal tip
            if (child.isLeaf() && branchLength != 0.0) {
                if (popSE) {
                    nodeMath.setVarianceForTip(childIdx,nodeMath.getVarianceForNode(childIdx) + 1);
                }

                PruneLikelihoodUtils.populateACEf(nodeMath, nodeMath.getVarianceForNode(childIdx), nTraits, childIdx);

                calculateLmrForTips(nodeMath, traitValuesArr, nTraits, childIdx);

                // add up to this node
                thisNodeL += nodeMath.getLForNode(childIdx);
                thisNodeR += nodeMath.getRForNode(childIdx);
                MatrixUtilsContra.vectorAdd(thisNodeMVec, nodeMath.getTempVec(), thisNodeMVec);

            } else {
                // if child is an sampled ancestor
                // it will not be included
                if (!child.isDirectAncestor()) {

                    // (2) child is a normal internal node, i.e. child does not have a sampled ancestor
                    if (!child.getChild(0).isDirectAncestor() && !child.getChild(1).isDirectAncestor()) {

                        PruneLikelihoodUtils.populateACEf(nodeMath, branchLength, nTraits, childIdx);

                        pruneNode(child, nTraits, traitValuesArr, pcmc, nodeMath, popSE);

                        calculateLmrForInternalNodes(nodeMath, nTraits, childIdx);

                        // add up to this node
                        thisNodeR += nodeMath.getRForNode(childIdx);
                        MatrixUtilsContra.vectorAdd(thisNodeMVec, nodeMath.getTempVec(), thisNodeMVec);
                        thisNodeL += nodeMath.getLForNode(childIdx);

                    } else {
                        // (3) child is an internal node and has a sampled ancestor below
                        // child will be used as a normal tip
                        // trait values of the sampled ancestor will used to calculate AbCdEf and Lmr
                        Node gcSA = child.getChild(0);
                        if(child.getChild(1).isDirectAncestor()) {
                            gcSA = child.getChild(1);
                        }
                        int gcSANr = gcSA.getNr();

                        PruneLikelihoodUtils.populateACEf(nodeMath, branchLength, nTraits, gcSANr);

                        calculateLmrForTips(nodeMath, traitValuesArr, nTraits, gcSANr);

                        // add up to this node
                        thisNodeL += nodeMath.getLForNode(gcSANr);
                        thisNodeR += nodeMath.getRForNode(gcSANr);
                        MatrixUtilsContra.vectorAdd(thisNodeMVec, nodeMath.getTempVec(), thisNodeMVec);

                        // prune the subtree rooted at the child node
                        pruneNode(child, nTraits, traitValuesArr, pcmc, nodeMath, popSE);
                        if (nodeMath.isSingularMatrix()) {
                            nodeMath.setLikelihoodForSampledAncestors(Double.NEGATIVE_INFINITY);
                        } else {
                            // get the trait values the child node, i.e. the trait values at the sampled ancestor
                            nodeMath.setTraitsVecForSampledAncestor(traitValuesArr, gcSANr);
                            
                            // calculate the likelihood rooted at this sampled ancestor
                            double logPSA; // currently, we use this flag to tell if shrinkage method is used to distinguish the calculation of likelihood
                            if(popSE) {
                                logPSA =
                                        MatrixUtilsContra.vecTransScalarMultiply(
                                                nodeMath.getSampledAncestorTraitsVec(),
                                                nodeMath.getLForNode(childIdx),
                                                nTraits) +
                                                MatrixUtilsContra.vectorDotMultiply(
                                                        nodeMath.getSampledAncestorTraitsVec(),
                                                        nodeMath.getMVecForNode(childIdx)) +
                                                nodeMath.getRForNode(childIdx);
                            } else {
                                logPSA = nodeMath.getLForNode(childIdx) *
                                        MatrixUtilsContra.tVecDotMatrixDotVec(
                                                nodeMath.getSampledAncestorTraitsVec(),
                                                nodeMath.getTraitRateMatrixInverse(),
                                                nTraits) +
                                        MatrixUtilsContra.vectorDotMultiply(
                                                nodeMath.getSampledAncestorTraitsVec(),
                                                nodeMath.getMVecForNode(childIdx)) +
                                        nodeMath.getRForNode(childIdx);
                            }
                            // add up to the likelihood for sampled ancestors in the tree
                            nodeMath.setLikelihoodForSampledAncestors(
                                    nodeMath.getLikelihoodForSampledAncestors() +
                                            logPSA);
                        }
                    }
                }
            }
        }
        if(nodeMath.hasMissingDataSpecies(thisNodeIdx)) {
            int speciesHadDataIndex = nodeMath.getSpeciesToIgnoreIndex(thisNodeIdx);
            double compoundVariance = nodeMath.getVarianceForNode(thisNodeIdx) + nodeMath.getVarianceForNode(speciesHadDataIndex);
            nodeMath.setVarianceForTip(thisNodeIdx, compoundVariance);
            nodeMath.setExpectationForIntNode(thisNodeIdx, nodeMath.getExpectationForNode(speciesHadDataIndex));
        } else {
            // estimate the node variance
            nodeMath.setVarianceForParent(thisNodeIdx, node.getLength() * pcmc.getRateForBranch(node), node.getChild(0).getNr(), node.getChild(1).getNr());

            // estimate the node expectations
            nodeMath.setExpectationForParent(thisNodeIdx, node.getChild(0).getNr(), node.getChild(1).getNr());
        }
        // set lMat, mVec and r for this node
        nodeMath.setLForNode(thisNodeIdx, thisNodeL);
        nodeMath.setMVecForNode(thisNodeIdx, thisNodeMVec);
        nodeMath.setRForNode(thisNodeIdx, thisNodeR);
    }

    protected void calculateLmrForTips(NodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {}

    protected void calculateLmrForInternalNodes(NodeMath nodeMath, int nTraits, int nodeIdx) {}

    protected double calculateLikelihood(NodeMath nodeMath, double l0, double[] m0, double r0, int rootIdx) {return 1.0;}

    /*
     * This function populates the arrays for species having missing data
     */
    private void prune(Node node){
        List<Node> children = node.getChildren();
        for (Node child : children) {
            if(!child.isLeaf()){
                prune(child);
                Node c1 = child.getChild(0);
                Node c2 = child.getChild(1);
                // (1) an internal node having both tips to ignore
                // -> this internal node will be ignored
                // -> its parent will be set as a node having missing data
                // -> indicating that its parent will get its sibling during likelihood calculation
                if(nodeMath.isSpeciesToIgnore(c1.getNr()) && nodeMath.isSpeciesToIgnore(c2.getNr())){
                    nodeMath.setSpeciesToIgnore(child.getNr());
                    nodeMath.setNodeHasMissingData(node.getNr());

                    Node sib = node.getChild(0);
                    int sibNr = sib.getNr();
                    if(sibNr  == child.getNr()){
                        nodeMath.setSpeciesToIgnoreIndex(node.getNr(), node.getChild(1).getNr());
                    } else {
                        nodeMath.setSpeciesToIgnoreIndex(node.getNr(), sibNr);
                    }
                }
                // (2) one of the child of the internal node is a species to ignore
                // -> set this internal node as a node having missing data
                // -> this internal node will get the child having data during likelihood calculation
                if(nodeMath.isSpeciesToIgnore(c1.getNr()) && !nodeMath.isSpeciesToIgnore(c2.getNr())){
                    nodeMath.setNodeHasMissingData(child.getNr());
                    nodeMath.setSpeciesToIgnoreIndex(child.getNr(), c2.getNr());
                }

                if(!nodeMath.isSpeciesToIgnore(c1.getNr()) && nodeMath.isSpeciesToIgnore(c2.getNr())){
                    nodeMath.setNodeHasMissingData(child.getNr());
                    nodeMath.setSpeciesToIgnoreIndex(child.getNr(), c1.getNr());
                }


            }
        }

    }
    @Override
    public List<String> getArguments() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public List<String> getConditions() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        // TODO Auto-generated method stub

    }
}
