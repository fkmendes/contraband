package contraband.prunelikelihood;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.*;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import contraband.valuewrappers.OneValueContTraits;
import org.apache.commons.math3.linear.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

@Description("This class implements likelihood for continuous traits under Ornstein–Uhlenbeck process.\n" +
        "The calculation uses Venelin's PCM likelihood.")


public class JOUPruneLikelihood extends OUPruneLikelihoodProcess {


    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    @Override
    public double calculateLogP() {

        return logP;
    }

    @Override
    protected RealMatrix calculatePhiMatrix (Node node, OUNodeMath nodeMath) {
        return OUPruneUtils.getPhiRM(node, nodeMath.getAlphaMatrix());
    }

    @Override
    protected RealVector calculateOmegaVector (Node node, OUNodeMath nodeMath, RealMatrix phiRM) {
        int xi = nodeMath.getJumpForNode(node.getNr());

        RealVector iMinusPhiTheta = OUPruneUtils.getOmegaVec(nodeMath.getThetaForNode(node.getNr()), phiRM, nodeMath.getIdentityMatrix());

        if(xi == 0) {
            return iMinusPhiTheta;
        } else{
            return OUPruneUtils.getOmegaVec(nodeMath.getJumpMeanVec(), phiRM, iMinusPhiTheta);
        }
    }

}
