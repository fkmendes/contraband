package contraband.prunelikelihood;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;

import java.util.List;

public class OrderedDiscreteTraits extends ThresholdModel {
    final public Input<RealParameter> thresholdsInput = new Input<>("threshold", "Thresholds for mapping ordered discrete traits.", Input.Validate.REQUIRED);

    private int nSpecies;
    private int nTraits;

    @Override
    public void initAndValidate() {

        super.initAndValidate();
        nSpecies = getSpeciesNr();
        nTraits = getTraitNr();
        if (liabilitiesInput.get().getDimension() != nSpecies * nTraits){
            liabilitiesInput.get().setDimension(nSpecies * nTraits);
        }


    }

}
