package contraband;

import java.util.List;
import java.util.Random;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;

public class RateCategoryIdentifiability extends Distribution {

	final public Input<RealParameter> rateValuesInput = new Input<>("rates", "the rate parameters associated with each category.", Input.Validate.REQUIRED);
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
	}
	
	@Override
	public double calculateLogP() {	
		Double[] rates = rateValuesInput.get().getValues();
		
		boolean ratesAreGo = true; // can be optima
		
		double lastRateValue = Double.NEGATIVE_INFINITY;
		for (double thetaValue: rates) {
			if (thetaValue < lastRateValue) {
				ratesAreGo = false;
			} else {
				lastRateValue = thetaValue;
			}
		}	
		
		if (ratesAreGo) { return 0.0; }
		else { return Double.NEGATIVE_INFINITY; }
	}

	@Override
	protected boolean requiresRecalculation() {
		return super.requiresRecalculation();
	}

	@Override
	public List<String> getArguments() {
		return null;
	}

	@Override
	public List<String> getConditions() {
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		;
	}
	
}

