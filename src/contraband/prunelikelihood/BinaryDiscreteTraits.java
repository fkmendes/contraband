package contraband.prunelikelihood;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;

import java.util.List;

public class BinaryDiscreteTraits extends ThresholdModel{

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    @Override
    public double calculateLogP() {

        return logP;
    }


}
