package contraband.prunelikelihood;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Tree;
import beast.evolution.tree.Node;
import contraband.clock.RateCategoryClockModel;
import contraband.math.MatrixUtilsContra;
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
    private OUNodeMath nodeMath;


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

        nodeMath.populateAlphaMatrix();
        nodeMath.performAlphaDecompostion();

        nodeMath.setSingularMatrix(false);

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
        logP = l0Mat.preMultiply(rootValuesVec).dotProduct(rootValuesVec) + rootValuesVec.dotProduct(m0Vec) + r0;

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

                // For OU, variance matrix, Phi and Omega need to be calculated for this node.
                RealMatrix phiRM = calculatePhiMatrix(child, nodeMath);
                RealVector omegaVec = calculateOmegaVector(child, nodeMath, phiRM);

                // inverse of variance-covariance of this node
                nodeMath.populateVarianceCovarianceMatrix(child);

                // For OU
                // matrix A, C, E
                // vector b, d
                // and value f
                // need to be calculated
                // for this node
                RealMatrix aMat = OUPruneUtils.getAMatForOU(nodeMath.getInverseVarianceMatrix());
                RealMatrix eMat = OUPruneUtils.getEMatForOU(phiRM, nodeMath.getInverseVarianceMatrix());
                RealMatrix cMat = OUPruneUtils.getCMatForOU(phiRM, eMat);
                double f = OUPruneUtils.getFforOU(omegaVec, nodeMath.getInverseVarianceMatrix(), nodeMath.getVCVMatDet(), nTraits);
                RealVector bVec = OUPruneUtils.getBVecForOU(nodeMath.getInverseVarianceMatrix(), omegaVec);
                RealVector dVec = OUPruneUtils.getDVecForOU(eMat, omegaVec);
                nodeMath.setAbCdEfOmegaPhiForNode(childIdx, aMat, bVec, cMat, dVec, eMat, f, omegaVec, phiRM);

                if (child.isLeaf()) {
                    // vector of trait values at this tip
                    RealVector traitsVec = traitValuesList.get(childIdx);

                    // add up L matrix
                    nodeMath.setLMatForNode(childIdx, OUPruneUtils.getLMatForOULeaf(cMat));
                    thisNodeLMat = thisNodeLMat.add(OUPruneUtils.getLMatForOULeaf(cMat));

                    // add up r value
                    nodeMath.setRForNode(childIdx, OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f));
                    thisNodeR += OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f);

                    // add up m vector
                    nodeMath.setMVecForNode(childIdx, OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec));
                    thisNodeMVec = thisNodeMVec.add(OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec));
                } else {

                    pruneOU(child, nodeMath);

                    nodeMath.performAPlusOperations(child, aMat);
                    // (aMat + lMat).inverse
                    RealMatrix aPlusLInv = nodeMath.getAPlusLInv();
                    // determinant of -2 * (aMat + lMat) in log space
                    double logDetVNode = nodeMath.getNegativeTwoAPlusLDet();

                    // add up r value
                    double r = OUPruneUtils.getRForOUIntNode(bVec, nodeMath.getMVecForNode(childIdx), aPlusLInv, f, nodeMath.getRForNode(childIdx), nTraits, logDetVNode);
                    nodeMath.setRForNode(childIdx, r);
                    thisNodeR += r;

                    RealMatrix eAPlusLInv = eMat.multiply(aPlusLInv);

                    // add up m vector
                    RealVector m = OUPruneUtils.getMVecForOUIntNode(eAPlusLInv, bVec, nodeMath.getMVecForNode(childIdx), dVec);
                    nodeMath.setMVecForNode(childIdx, m);
                    thisNodeMVec = thisNodeMVec.add(m);

                    // add up L matrix
                    RealMatrix l = OUPruneUtils.getLMatForOUIntNode(cMat, eMat, eAPlusLInv);
                    nodeMath.setLMatForNode(childIdx, l);
                    thisNodeLMat = thisNodeLMat.add(l);
                }
            }
            // set L, m , r for this node
            nodeMath.setLMatForNode(thisNodeIdx, thisNodeLMat);
            nodeMath.setMVecForNode(thisNodeIdx, thisNodeMVec);
            nodeMath.setRForNode(thisNodeIdx, thisNodeR);
    }

    protected RealMatrix calculatePhiMatrix (Node node, OUNodeMath nodeMath) { return null;}

    protected RealVector calculateOmegaVector (Node node, OUNodeMath nodeMath, RealMatrix PhiMat) { return null;}

    private void populateTraitValuesList(KeyRealParameter traitValues, Tree tree, List<RealVector> traitValuesList) {
        // according to node number of tips
        for (int i = 0; i < tree.getLeafNodeCount(); i ++) {
            // get all traits values for this species
            Double[] traitForSpecies = traitValues.getRowValues(tree.getNode(i).getID());
            traitValuesList.add(i, new ArrayRealVector(traitForSpecies));
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
