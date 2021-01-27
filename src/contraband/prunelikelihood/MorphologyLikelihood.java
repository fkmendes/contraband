package contraband.prunelikelihood;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import contraband.math.GeneralNodeMath;
import contraband.math.MatrixUtilsContra;
import contraband.utils.MorphologyLikelihoodUtils;
import java.util.List;
import java.util.Random;

public class MorphologyLikelihood extends Distribution {
    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Input.Validate.REQUIRED);
    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel", "the rate or optimum on each branch", Input.Validate.REQUIRED);
    final public Input<GeneralNodeMath> nodeMathInput = new Input<>("nodeMath","Node information that will be used in PCM likelihood calculation.", Input.Validate.REQUIRED);
    final public Input<MorphologicalData> traitInput = new Input<>("trait","Morphological data set.", Input.Validate.REQUIRED);

    private Tree tree;
    private BranchRateModel.Base branchRateModel;
    private double[] traitValuesArr;
    private GeneralNodeMath nodeMath;

    private int nTraits;
    private boolean transformData;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // get the tree
        tree = treeInput.get();

        // get clock model
        branchRateModel = branchRateModelInput.get();

        nodeMath = nodeMathInput.get();

        transformData = traitInput.get().getTransformDataFlag();

        int nSpecies = traitInput.get().getSpeciesNr();
        nTraits = traitInput.get().getTotalTraitNr();
        traitValuesArr = new double[nSpecies * nTraits];
        traitValuesArr = traitInput.get().getMorphologicalData();
    }


    @Override
    public double calculateLogP() {
        // reject the tree with a sampled ancestor being the child of root
        if(tree.getRoot().getChild(0).isDirectAncestor() || tree.getRoot().getChild(1).isDirectAncestor()) {
            return Double.NEGATIVE_INFINITY;
        }

        // initiate the likelihood for sampled ancestors
        nodeMath.setLikelihoodForSampledAncestors(0.0);

        // prune the tree by starting from the root
        prune(tree.getRoot(), nTraits, traitValuesArr, branchRateModel, nodeMath);


        // if at some internal node, (aMat + lMat) is singular or -2 * (aMat + lMat) is singular,
        // when calculating (aMat + lMat).inverse and det[-2 * (aMat + lMat)],
        // reject the state
        if (nodeMath.getSingularMatrixFlag()) {
            return Double.NEGATIVE_INFINITY;
        }

        int rootIdx = tree.getRoot().getNr();
        // populate the root values
        nodeMath.populateRootValuesVec(rootIdx);

        // get L, m and r at the root
        double l0 = nodeMath.getLForNode(rootIdx);
        double[] m0 = nodeMath.getMVecForNode(rootIdx);
        double r0 = nodeMath.getRForNode(rootIdx);

        // calculate likelihood in log space
        logP = calculateLikelihood(nodeMath, l0, m0, r0, nTraits, rootIdx);

        return logP;
    }

    public void prune(Node node, int nTraits, double[] traitValuesArr,
                      BranchRateModel.Base pcmc, GeneralNodeMath nodeMath) {

        int thisNodeIdx = node.getNr();

        // initialize before sum
        double thisNodeL = 0.0;
        double [] thisNodeMVec = nodeMath.getInitMVec().clone();
        double thisNodeR = 0.0;

        List<Node> children = node.getChildren();

        for (Node child : children) {

            int childIdx = child.getNr();

            double branchLength = child.getLength() * pcmc.getRateForBranch(child);
            //nodeMath.setVarianceForTip(childIdx, branchLength);

            // (1) child is a normal tip
            if (child.isLeaf() && branchLength != 0.0) {

                //nodeMath.setVarianceForTip(childIdx,nodeMath.getVarianceForNode(childIdx) + 1);

                //PruneLikelihoodUtils.populateACEf(nodeMath, nodeMath.getVarianceForNode(childIdx), nTraits, childIdx);
                populateAbCdEfForNode(nodeMath, branchLength, nTraits, childIdx);

                populateLmrForTips(nodeMath, traitValuesArr, nTraits, childIdx);

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

                        //PruneLikelihoodUtils.populateACEf(nodeMath, branchLength, nTraits, childIdx);
                        populateAbCdEfForNode(nodeMath, branchLength, nTraits, childIdx);

                        prune(child, nTraits, traitValuesArr, pcmc, nodeMath);

                        populateLmrForInternalNodes(nodeMath, nTraits, childIdx);

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

                        //PruneLikelihoodUtils.populateACEf(nodeMath, branchLength, nTraits, gcSANr);
                        populateAbCdEfForNode(nodeMath, branchLength, nTraits, gcSANr);

                        populateLmrForTips(nodeMath, traitValuesArr, nTraits, gcSANr);

                        // add up to this node
                        thisNodeL += nodeMath.getLForNode(gcSANr);
                        thisNodeR += nodeMath.getRForNode(gcSANr);
                        MatrixUtilsContra.vectorAdd(thisNodeMVec, nodeMath.getTempVec(), thisNodeMVec);

                        // prune the subtree rooted at the child node
                        prune(child, nTraits, traitValuesArr, pcmc, nodeMath);
                        if (nodeMath.getSingularMatrixFlag()) {
                            nodeMath.setLikelihoodForSampledAncestors(Double.NEGATIVE_INFINITY);
                        } else {
                            // get the trait values the child node, i.e. the trait values at the sampled ancestor
                            nodeMath.setTraitsVecForSampledAncestor(traitValuesArr, gcSANr);

                            // calculate the likelihood rooted at this sampled ancestor
                            double logPSA = calculateLikelihoodForSA(nodeMath, childIdx);

                            // add up to the likelihood for sampled ancestors in the tree
                            nodeMath.setLikelihoodForSampledAncestors(
                                    nodeMath.getLikelihoodForSampledAncestors() +
                                            logPSA);
                        }
                    }
                }
            }
        }

        // estimate the node variance
        //nodeMath.setVarianceForParent(thisNodeIdx, node.getLength() * pcmc.getRateForBranch(node), node.getChild(0).getNr(), node.getChild(1).getNr());

        // estimate the node expectations
        //nodeMath.setExpectationForParent(thisNodeIdx, node.getChild(0).getNr(), node.getChild(1).getNr());

        // set lMat, mVec and r for this node
        nodeMath.setLForNode(thisNodeIdx, thisNodeL);
        nodeMath.setMVecForNode(thisNodeIdx, thisNodeMVec);
        nodeMath.setRForNode(thisNodeIdx, thisNodeR);
    }

    private void populateAbCdEfForNode (GeneralNodeMath nodeMath, double branchLength, int nTraits, int nodeIdx) {
            MorphologyLikelihoodUtils.populateACEf(nodeMath, branchLength, nTraits, nodeIdx);
    }

    private void populateLmrForTips(GeneralNodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        if(transformData) {
            MorphologyLikelihoodUtils.populateLmrForTipTransform(nodeMath, traitValuesArr, nTraits, nodeIdx);
        } else {
            MorphologyLikelihoodUtils.populateLmrForTip(nodeMath, traitValuesArr, nTraits, nodeIdx);
        }
    }


    private void populateLmrForInternalNodes(GeneralNodeMath nodeMath, int nTraits, int nodeIdx) {
        if(transformData) {
            MorphologyLikelihoodUtils.populateLmrForInternalNodeTransform(nodeMath, nTraits, nodeIdx);
        } else {
            MorphologyLikelihoodUtils.populateLmrForInternalNode(nodeMath, nTraits, nodeIdx);
        }
    }

    private double calculateLikelihood(GeneralNodeMath nodeMath, double l0, double[] m0, double r0, int nTraits, int rootIdx){
        if(transformData){
            return MatrixUtilsContra.vecTransScalarMultiply(nodeMath.getRootValuesArr(),
                    l0, nTraits) +
                    MatrixUtilsContra.vectorDotMultiply(nodeMath.getRootValuesArr(), m0) +
                    r0 +
                    nodeMath.getLikelihoodForSampledAncestors();
        } else {
            return l0 * MatrixUtilsContra.tVecDotMatrixDotVec(
                    nodeMath.getRootValuesArr(),
                    nodeMath.getTraitRateMatrixInverse(),
                    nTraits) +
                    MatrixUtilsContra.vectorDotMultiply(
                            nodeMath.getRootValuesArr(),
                            m0) +
                    r0 + nodeMath.getLikelihoodForSampledAncestors();
        }
    }

    private double calculateLikelihoodForSA (GeneralNodeMath nodeMath, int childIdx) {
        if(transformData) {
            return MatrixUtilsContra.vecTransScalarMultiply(
                    nodeMath.getSampledAncestorTraitsVec(),
                    nodeMath.getLForNode(childIdx),
                    nTraits) +
                    MatrixUtilsContra.vectorDotMultiply(
                            nodeMath.getSampledAncestorTraitsVec(),
                            nodeMath.getMVecForNode(childIdx)) +
                    nodeMath.getRForNode(childIdx);
        } else {
            return nodeMath.getLForNode(childIdx) *
                    MatrixUtilsContra.tVecDotMatrixDotVec(
                            nodeMath.getSampledAncestorTraitsVec(),
                            nodeMath.getTraitRateMatrixInverse(),
                            nTraits) +
                    MatrixUtilsContra.vectorDotMultiply(
                            nodeMath.getSampledAncestorTraitsVec(),
                            nodeMath.getMVecForNode(childIdx)) +
                    nodeMath.getRForNode(childIdx);
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
