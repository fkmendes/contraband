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
        nodeMath.populateAlphaMatrix();
        nodeMath.performAlphaDecomposition(nodeMath.getAlphaMatrix());

        super.populateLogP();

        return getLogP();
    }
}
