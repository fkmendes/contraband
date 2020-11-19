package contraband.prunelikelihood;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.*;
import beast.evolution.tree.Tree;
import beast.evolution.tree.Node;
import contraband.math.MatrixUtilsContra;
import contraband.utils.PruneLikelihoodUtils;
import contraband.valuewrappers.OneValueContTraits;
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
    private double[] traitValuesArr;
    private OUNodeMath nodeMath;


    @Override
    public void initAndValidate() {
        //super.initAndValidate();
        tree = treeInput.get();

        KeyRealParameter traitsValues = traitsValuesInput.get();

        nTraits = traitsValues.getMinorDimension1();
        int nSpecies = traitsValues.getMinorDimension2();

        traitValuesArr = new double[nSpecies * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitsValues, tree, nTraits, traitValuesArr);

        nodeMath = nodeMathInput.get();
    }

    protected void populateLogP() {

        nodeMath.populateAlphaMatrix();

        nodeMath.performAlphaDecompostion();


        pruneOU(tree.getRoot(), nodeMath);

        RealMatrix l0Mat = nodeMath.getLMat();
        RealVector m0Vec = nodeMath.getMVec();
        double r0 = nodeMath.getR();

        RealVector rootValuesVec = nodeMath.getRootValuesVec();

        // calculate likelihood
        // loglik <- X0 * L0 * X0+ m0 * X0 + r0
        logP = l0Mat.preMultiply(rootValuesVec).dotProduct(rootValuesVec) + rootValuesVec.dotProduct(m0Vec) + r0;
    }

    // getters
    public double getLogP () { return logP; }


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
                //RealMatrix phiRM = getPhiRM(child, alphaMat);
                RealMatrix phiRM = calculatePhiMatrix(child, nodeMath);

                //RealVector omegaVec = getOmegaVec(thetaVec, phiRM, identity);
                RealVector omegaVec = calculateOmegaVector(nodeMath, phiRM);

                // inverse of variance-covariance of this node
                /*
                RealMatrix invVCVMat = getInverseVarianceRMForOU(child, vcvMatDetArr,
                        sigmaRM, sigmaeRM,
                        pMat, inverseP, decompositionH, nTraits);
                 double varianceRMDet = vcvMatDetArr[childIdx];
                 */
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
                double f = OUPruneUtils.getFforOU(omegaVec, nodeMath.getInverseVarianceMatrix(), nodeMath.getDetVarianceMatrix(), nTraits);
                RealVector bVec = OUPruneUtils.getBVecForOU(nodeMath.getInverseVarianceMatrix(), omegaVec);
                RealVector dVec = OUPruneUtils.getDVecForOU(eMat, omegaVec);

                if (child.isLeaf()) {
                    // vector of trait values at this tip
                    //RealVector traitsVec = traitsValuesList.get(childIdx);
                    double[] traitsArr = new double [nTraits];
                    MatrixUtilsContra.getMatrixRow(traitValuesArr, childIdx, nTraits, traitsArr);
                    RealVector traitsVec = new ArrayRealVector(traitsArr);

                    // set the L matrix
                    thisNodeLMat = thisNodeLMat.add(OUPruneUtils.getLMatForOULeaf(cMat));

                    // set r value
                    thisNodeR += OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f);

                    // set m vector
                    thisNodeMVec = thisNodeMVec.add(OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec));
                } else {

                    pruneOU(child, nodeMath);

                    // (aMat + lMat).inverse
                    //RealMatrix aPlusLInv = OUPruneUtils.getInvAPlusLRM(child, negativeTwoAplusLDetArr, aMat, lMatList.get(childIdx));

                    // determinant of -2 * (aMat + lMat)
                    //double logDetVNode = OUPruneUtils.negativeTwoAplusLDetArr[childIdx];

                    // set r value
                    //thisNodeR += OUPruneUtils.getRForOUIntNode(bVec, mVecList.get(childIdx), aPlusLInv, f, rArr[childIdx], nTraits, logDetVNode);

                    //RealMatrix eAPlusLInv = eMat.multiply(aPlusLInv);
                    // set m vector
                    //thisNodeMVec = thisNodeMVec.add(OUPruneUtils.getMVecForOUIntNode(eAPlusLInv, bVec, mVecList.get(childIdx), dVec));

                    // set L matrix
                    //thisNodeLMat = thisNodeLMat.add(OUPruneUtils.getLMatForOUIntNode(cMat, eMat, eAPlusLInv));
                }
            }
            //lMatList.set(thisNodeIdx, thisNodeLMat);
            //mVecList.set(thisNodeIdx, thisNodeMVec);
            //rArr[thisNodeIdx] = thisNodeR;
    }

    protected RealMatrix calculatePhiMatrix (Node node, OUNodeMath nodeMath) { return null;}

    protected RealVector calculateOmegaVector (OUNodeMath nodeMath, RealMatrix PhiMat) { return null;}

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
