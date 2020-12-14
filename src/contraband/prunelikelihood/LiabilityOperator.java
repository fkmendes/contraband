package contraband.prunelikelihood;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import outercore.parameter.KeyRealParameter;

public class LiabilityOperator extends Operator {
    final public Input<RealParameter> liabilitiesInput = new Input<>("liability", "Continuous random variables underneath observed traits.", Input.Validate.REQUIRED);
    final public Input<KeyRealParameter> traitsValuesInput = new Input<>("traits","Trait values at tips.", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {

    }

    @Override
    public double proposal() {

        return Double.POSITIVE_INFINITY;
    }
}
