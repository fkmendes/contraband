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


public class DOUPruneLikelihood extends OUPruneLikelihoodProcess {


    @Override
    public void initAndValidate() {
        super.initAndValidate();
        nodeMath.initValidateDOU();
    }

    @Override
    public double calculateLogP() {
        nodeMath.populateAlphaMatrix();
        nodeMath.updateDOUParameters();
        nodeMath.performAlphaDecomposition(nodeMath.getDAlphaMatrix());

        super.populateLogP();

        return getLogP();
    }
}
