package contraband.prunelikelihood;

import beast.core.Description;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import org.apache.commons.math3.linear.*;


@Description("This class implements likelihood for continuous traits under Ornstein–Uhlenbeck process.\n" +
        "The calculation uses Venelin's PCM likelihood.")

public class OUPruneLikelihood extends OUPruneLikelihoodProcess {

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    @Override
    public double calculateLogP() {
        super.populateLogP();

        return getLogP();
    }

    @Override
    protected RealMatrix calculatePhiMatrix (Node node, OUNodeMath nodeMath) {
        return OUPruneUtils.getPhiRM(node, nodeMath.getAlphaMatrix());
    }

    @Override
    protected RealVector calculateOmegaVector (Node node, OUNodeMath nodeMath, RealMatrix phiRM) {
        return OUPruneUtils.getOmegaVec(nodeMath.getThetaForNode(node.getNr()), phiRM, nodeMath.getIdentityMatrix());
    }


}
