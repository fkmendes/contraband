package contraband.otherlikelihood;

import java.util.List;
import java.util.Random;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;

/*
 * This is a "hack" likelihood that helps us reject colors (e.g., rates, adaptive optima) that are not in increasing order
 */
public class RateCategoryIdentifiability extends Distribution {

	final public Input<RealParameter> rateValuesInput = new Input<>("rates", "the rate parameters associated with each category.", Input.Validate.REQUIRED);
	
	Double[] rates;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
	}
	
	@Override
	public double calculateLogP() {	
		rates = rateValuesInput.get().getValues();
		
		boolean ratesAreGo = true; // can be optima
		
		double lastRateValue = Double.NEGATIVE_INFINITY;
		for (double thetaValue: rates) {
			if (thetaValue < lastRateValue) {
				ratesAreGo = false;
			} else {
				lastRateValue = thetaValue;
			}
		}	
		
		if (ratesAreGo) { logP = 0.0; }
		else { logP = Double.NEGATIVE_INFINITY; }
		return logP;
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

