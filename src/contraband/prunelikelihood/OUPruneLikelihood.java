package contraband.prunelikelihood;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.*;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import contraband.utils.PruneLikelihoodUtils;
import contraband.valuewrappers.OneValueContTraits;
import org.apache.commons.math3.linear.*;
import outercore.parameter.KeyRealParameter;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

@Description("This class implements likelihood for continuous traits under Ornstein–Uhlenbeck process.\n" +
        "The calculation uses Venelin's PCM likelihood.")


public class OUPruneLikelihood extends OUPruneLikelihoodProcess {

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    @Override
    protected RealMatrix calculatePhiMatrix (Node node, OUNodeMath nodeMath) {
        return OUPruneUtils.getPhiRM(node, nodeMath.getAlphaMatrix());
    }

    @Override
    protected RealVector calculateOmegaVector (OUNodeMath nodeMath, RealMatrix phiRM) {
        return OUPruneUtils.getOmegaVec(nodeMath.getThetaVec(), phiRM, nodeMath.getIdentityMatrix());
    }

    @Override
    public double calculateLogP() {
        super.populateLogP();
        return getLogP();
    }
}
