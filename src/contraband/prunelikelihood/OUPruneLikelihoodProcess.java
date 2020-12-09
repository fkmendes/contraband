package contraband.prunelikelihood;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.evolution.tree.Tree;
import beast.evolution.tree.Node;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.utils.PruneLikelihoodUtils;
import org.apache.commons.math3.linear.*;
import outercore.parameter.KeyRealParameter;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

@Description("This class implements likelihood for continuous traits under Ornstein–Uhlenbeck process.\n" +
        "The calculation uses Venelin's PCM likelihood.")

public abstract class OUPruneLikelihoodProcess extends Distribution {
    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Input.Validate.REQUIRED);
    final public Input<KeyRealParameter> traitsValuesInput = new Input<>("traits","Trait values at tips.", Input.Validate.REQUIRED);
    final public Input<OUNodeMath> nodeMathInput = new Input<>("nodeMath","Parameters at each node.", Input.Validate.REQUIRED);

    private Tree tree;
    private int nTraits;
    private List<RealVector>  traitValuesList;
    protected OUNodeMath nodeMath;


    @Override
    public void initAndValidate() {
        super.initAndValidate();
        tree = treeInput.get();

        KeyRealParameter traitsValues = traitsValuesInput.get();

        nTraits = traitsValues.getMinorDimension1();
        int nSpecies = traitsValues.getMinorDimension2();

        traitValuesList = new ArrayList<>(nSpecies);
        populateTraitValuesList(traitsValues, tree, traitValuesList);

        nodeMath = nodeMathInput.get();
    }

    protected void populateLogP() {
        nodeMath.updateSigmaMatrix();
        nodeMath.updateThetaVectors();
        nodeMath.updateRootValues();

        nodeMath.setSingularMatrix(false);
        nodeMath.setLikelihoodForSA(0.0);

        pruneOU(tree.getRoot(), nodeMath);
    }

    // getters
    public double getLogP () {
        if (nodeMath.getSingularMatrix()) {
            return logP = Double.NEGATIVE_INFINITY;
        }


        int rootIdx = tree.getRoot().getNr();
        RealMatrix l0Mat = nodeMath.getLMatForNode(rootIdx);
        RealVector m0Vec = nodeMath.getMVecForNode(rootIdx);
        double r0 = nodeMath.getRForNode(rootIdx);

        RealVector rootValuesVec = nodeMath.getRootValuesVec();

        // calculate likelihood
        // loglik <- X0 * L0 * X0+ m0 * X0 + r0
        logP = l0Mat.preMultiply(rootValuesVec).dotProduct(rootValuesVec) + rootValuesVec.dotProduct(m0Vec) + r0 + nodeMath.getLikelihoodForSA();

        return logP;
    }


    public void pruneOU (Node node, OUNodeMath nodeMath) {
        int thisNodeIdx = node.getNr();

        // initialize before sum
        RealMatrix thisNodeLMat = nodeMath.getIniMatrix().copy();
        RealVector thisNodeMVec = nodeMath.getIniVec().copy();
        double thisNodeR = 0.0;

            List<Node> children = node.getChildren();

            for (Node child : children) {

                int childIdx = child.getNr();

                if (child.isLeaf() && child.getLength() != 0.0) {
                    // vector of trait values at this tip
                    RealVector traitsVec = traitValuesList.get(childIdx);

                    // For OU, matrix A, C, E, vector b, d and value f need to be calculated for this node
                    populateAbCdEf(nodeMath, child, childIdx);

                    calculateLmrForTipNodes(nodeMath, traitsVec, childIdx);

                    // add up L matrix
                    thisNodeLMat = thisNodeLMat.add(nodeMath.getLMatForNode(childIdx));
                    // add up r value
                    thisNodeR += nodeMath.getRForNode(childIdx);
                    // add up m vector
                    thisNodeMVec = thisNodeMVec.add(nodeMath.getMVecForNode(childIdx));

                } else {
                    // if child is an sampled ancestor
                    // it will not be included
                    if (!child.isDirectAncestor()) {

                        // (2) child is a normal internal node, i.e. child does not have a sampled ancestor
                        if (!child.getChild(0).isDirectAncestor() && !child.getChild(1).isDirectAncestor()) {
                            // For OU, matrix A, C, E, vector b, d and value f need to be calculated for this node
                            populateAbCdEf(nodeMath, child, childIdx);

                            pruneOU(child, nodeMath);

                            calculateLmrForIntNode(child, nodeMath, childIdx);

                            // add up r value
                            thisNodeR += nodeMath.getRForNode(childIdx);
                            // add up m vector
                            thisNodeMVec = thisNodeMVec.add(nodeMath.getMVecForNode(childIdx));
                            // add up L matrix
                            thisNodeLMat = thisNodeLMat.add(nodeMath.getLMatForNode(childIdx));
                        }
                        else {
                            // (3) child is an internal node and has a sampled ancestor below
                            // child will be used as a normal tip
                            // trait values of the sampled ancestor will used to calculate AbCdEf and Lmr
                            Node gcSA = child.getChild(0);
                            if(child.getChild(1).isDirectAncestor()) {
                                gcSA = child.getChild(1);
                            }
                            int gcSANr = gcSA.getNr();

                            populateAbCdEf(nodeMath, gcSA, gcSANr);

                            // vector of trait values at this tip
                            RealVector saTraitsVec = traitValuesList.get(gcSANr);

                            calculateLmrForTipNodes(nodeMath, saTraitsVec, gcSANr);

                            // add up to this node
                            thisNodeLMat = thisNodeLMat.add(nodeMath.getLMatForNode(gcSANr));
                            thisNodeR += nodeMath.getRForNode(gcSANr);
                            thisNodeMVec = thisNodeMVec.add(nodeMath.getMVecForNode(gcSANr));

                            // prune the subtree rooted at the child node
                            pruneOU(child, nodeMath);
                            if (nodeMath.getSingularMatrix()) {
                                nodeMath.setLikelihoodForSA(Double.NEGATIVE_INFINITY);
                            } else {
                                // calculate the likelihood rooted at this sampled ancestor
                                double logPSA = nodeMath.getLMatForNode(childIdx).preMultiply(saTraitsVec).dotProduct(saTraitsVec) + saTraitsVec.dotProduct(nodeMath.getMVecForNode(childIdx)) + nodeMath.getRForNode(childIdx);

                                // add up to the likelihood for sampled ancestors in the tree
                                nodeMath.setLikelihoodForSA(
                                        nodeMath.getLikelihoodForSA() +
                                                logPSA);
                            }
                        }
                    }
                }
            }
            // set L, m , r for this node
            nodeMath.setLMatForNode(thisNodeIdx, thisNodeLMat);
            nodeMath.setMVecForNode(thisNodeIdx, thisNodeMVec);
            nodeMath.setRForNode(thisNodeIdx, thisNodeR);
    }

    protected RealMatrix calculatePhiMatrix (Node node, OUNodeMath nodeMath) {
        return OUPruneUtils.getPhiRM(node, nodeMath.getAlphaMatrix());
    }

    protected RealVector calculateOmegaVector (Node node, OUNodeMath nodeMath, RealMatrix phiRM) {
        return OUPruneUtils.getOmegaVec(nodeMath.getThetaForNode(node.getNr()), phiRM, nodeMath.getIdentityMatrix());
    }
    private void populateTraitValuesList(KeyRealParameter traitValues, Tree tree, List<RealVector> traitValuesList) {
        // according to node number of tips
        for (int i = 0; i < tree.getLeafNodeCount(); i ++) {
            // get all traits values for this species
            Double[] traitForSpecies = traitValues.getRowValues(tree.getNode(i).getID());
            traitValuesList.add(i, new ArrayRealVector(traitForSpecies));
        }
    }
    protected void populateAbCdEf (OUNodeMath nodeMath, Node node, int nodeIdx) {
        // For OU, variance matrix, Phi and Omega need to be calculated for this node.
        RealMatrix phiRM = calculatePhiMatrix(node, nodeMath);
        RealVector omegaVec = calculateOmegaVector(node, nodeMath, phiRM);

        // inverse of variance-covariance of this node
        nodeMath.populateVarianceCovarianceMatrix(node);
        nodeMath.setVarianceCovarianceMatrix(node, nodeIdx);

        RealMatrix aMat = OUPruneUtils.getAMatForOU(nodeMath.getInverseVarianceMatrix());
        RealMatrix eMat = OUPruneUtils.getEMatForOU(phiRM, nodeMath.getInverseVarianceMatrix());
        RealMatrix cMat = OUPruneUtils.getCMatForOU(phiRM, eMat);
        double f = OUPruneUtils.getFforOU(omegaVec, nodeMath.getInverseVarianceMatrix(), nodeMath.getVCVMatDet(), nTraits);
        RealVector bVec = OUPruneUtils.getBVecForOU(nodeMath.getInverseVarianceMatrix(), omegaVec);
        RealVector dVec = OUPruneUtils.getDVecForOU(eMat, omegaVec);

        nodeMath.setAbCdEfOmegaPhiForNode(nodeIdx, aMat, bVec, cMat, dVec, eMat, f, omegaVec, phiRM);
    }

    protected void calculateLmrForTipNodes(OUNodeMath nodeMath, RealVector traitsVec, int nodeIdx) {
        RealMatrix aMat = nodeMath.getAMatForNode(nodeIdx);
        RealMatrix cMat = nodeMath.getCMatForNode(nodeIdx);
        RealMatrix eMat = nodeMath.getEMatForNode(nodeIdx);
        RealVector bVec = nodeMath.getbVecForNode(nodeIdx);
        RealVector dVec = nodeMath.getdVecForNode(nodeIdx);
        double f = nodeMath.getfForNode(nodeIdx);

        // add up L matrix
        nodeMath.setLMatForNode(nodeIdx, OUPruneUtils.getLMatForOULeaf(cMat));
        //thisNodeLMat = thisNodeLMat.add(OUPruneUtils.getLMatForOULeaf(cMat));

        // add up r value
        nodeMath.setRForNode(nodeIdx, OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f));
        //thisNodeR += OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f);

        // add up m vector
        nodeMath.setMVecForNode(nodeIdx, OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec));
        //thisNodeMVec = thisNodeMVec.add(OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec));
    }

    protected void calculateLmrForIntNode(Node child, OUNodeMath nodeMath, int childIdx){
        nodeMath.performAPlusOperations(child, nodeMath.getAMatForNode(childIdx));
        // (aMat + lMat).inverse
        RealMatrix aPlusLInv = nodeMath.getAPlusLInv();
        // determinant of -2 * (aMat + lMat) in log space
        double logDetVNode = nodeMath.getNegativeTwoAPlusLDet();
        // eMat * (aMat + lMat).inverse
        RealMatrix eAPlusLInv = nodeMath.getEMatForNode(childIdx).multiply(aPlusLInv);

        double r = OUPruneUtils.getRForOUIntNode(nodeMath.getbVecForNode(childIdx), nodeMath.getMVecForNode(childIdx), aPlusLInv, nodeMath.getfForNode(childIdx), nodeMath.getRForNode(childIdx), nTraits, logDetVNode);
        nodeMath.setRForNode(childIdx, r);

        RealVector m = OUPruneUtils.getMVecForOUIntNode(eAPlusLInv, nodeMath.getbVecForNode(childIdx), nodeMath.getMVecForNode(childIdx), nodeMath.getdVecForNode(childIdx));
        nodeMath.setMVecForNode(childIdx, m);

        RealMatrix l = OUPruneUtils.getLMatForOUIntNode(nodeMath.getCMatForNode(childIdx), nodeMath.getEMatForNode(childIdx), eAPlusLInv);
        nodeMath.setLMatForNode(childIdx, l);
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
