package contraband.otherlikelihood;

import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import contraband.math.MVNUtils;
import contraband.valuewrappers.OneValueContTraits;

/*
 * This is the likelihood for the white noise (WN) model. This model assumes that the each species data comes from
 * a normal distribution. The normal distribution can be set to be the same for all species, or different for
 * all species ("Multiple regime" WN).
 * 
 * Note that under WN, we are basically saying that the trait data does not evolve along a tree (or, if we want to
 * think it in tree terms, it evolves along a star tree).
 */
public class WNLikelihood extends Distribution {

	final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	final public Input<RealParameter> sigmaSqsInput = new Input<>("sigmaSqs", "Sigma^2 of each species' normal density.", Validate.REQUIRED);
	final public Input<RealParameter> meansInput = new Input<>("mus", "Means of each species' normal density.", Validate.REQUIRED);
	final public Input<IntegerParameter> normalAssignmentsInput = new Input<>("normalAssignments", "Which normal density each species has.", Validate.REQUIRED);
	
	private OneValueContTraits sampleData;
	private String[] spNames;
	private Double[] mus, sigmaSqs, samples;
	
	// stored stuff
	private Double[] storedSamples;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		sampleData = oneTraitInput.get();
		spNames = sampleData.getSpNames();
		samples = sampleData.getTraitValues(0, spNames);
		
		mus = new Double[sampleData.getNSpp()];
		sigmaSqs = new Double[sampleData.getNSpp()];
		
		// stored stuff
		storedSamples = new Double[spNames.length]; 
	}
	
	@Override
	public double calculateLogP() {
		populateSampleData();
		
		Integer[] normalAssignments = normalAssignmentsInput.get().getValues();
		Double[] allSigmaSqsValues = sigmaSqsInput.get().getValues();
		Double[] allMusValues = meansInput.get().getValues();
		
		int i = 0;
		for (Integer assignment: normalAssignments) {
			mus[i] = allMusValues[assignment];
			sigmaSqs[i] = allSigmaSqsValues[assignment];
			i++;
		}
		
		logP = MVNUtils.getSampleMultipleNormalLogLk(samples, mus, sigmaSqs);
		
		return logP;
	}
	
	private void populateSampleData() {
		if (oneTraitInput.isDirty()) {
			sampleData = oneTraitInput.get();
			samples = sampleData.getTraitValues(0, spNames);
		}
	}
	
	@Override
	public boolean requiresRecalculation() {
		boolean dirty = true;
		return dirty;
	}
	
	@Override
	public void store() {
		System.arraycopy(samples, 0, storedSamples, 0, samples.length);
	}
	
	@Override
	public void restore() {		
		Double[] arrTmp;
		
		arrTmp = samples;
		samples = storedSamples;
		storedSamples = arrTmp;
	}
	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	
}

