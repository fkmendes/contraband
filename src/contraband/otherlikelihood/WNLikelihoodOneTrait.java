package contraband.otherlikelihood;

import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import contraband.utils.MVNUtils;

/**
 * @author Fabio K. Mendes
 */

@Description("This is the likelihood for the white noise (WN) model." +
		"This model assumes that the each species data comes from a normal distribution." +
		"The normal distribution can be set to be the same for all species, or different for all species." +
		"Note that under WN, we are basically saying that the trait data does not evolve along a tree" +
		"(or, if we want to think it in tree terms, it evolves along a star tree).")
public class WNLikelihoodOneTrait extends Distribution {

	final public Input<RealParameter> sigmaSqsInput = new Input<>("sigmaSqs", "Sigma^2 of each species' normal density.", Validate.REQUIRED);
	final public Input<RealParameter> meansInput = new Input<>("mus", "Means of each species' normal density.", Validate.REQUIRED);
	final public Input<IntegerParameter> normalAssignmentsInput = new Input<>("normalAssignments", "Which normal density each species has.", Validate.REQUIRED);
	final public Input<RealParameter> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED);
	// final public Input<OneValueContTraits> oneTraitInput = new Input<>("oneTraitData", "continuous data values for one trait.", Validate.REQUIRED); // original implementation (for the above line) with OneValueContTraits data wrapper
	
	// private OneValueContTraits sampleData; // original implementation with OneValueContTraits data wrapper
	private Double[] oneTraitDataArr;
	private Double[] mus, sigmaSqs;
	
	// stored stuff
	private Double[] storedOneTraitDataArr;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		RealParameter oneTraitValues = oneTraitInput.get();
		oneTraitDataArr = oneTraitValues.getValues();

		int nSpp = oneTraitValues.getMinorDimension2();
		mus = new Double[nSpp];
		sigmaSqs = new Double[nSpp];
		
		// stored stuff
		storedOneTraitDataArr = new Double[nSpp];
	}
	
	@Override
	public double calculateLogP() {

		// in case we ever sample tip values because of missing data, say
		boolean updateTipValues = false;
		if (oneTraitInput.isDirty()) { updateTipValues = true; }
		populateSampleData(updateTipValues);

		// TODO: later should add support for clocks and checking if operated on (if is dirty)
		Integer[] normalAssignments = normalAssignmentsInput.get().getValues();
		Double[] allSigmaSqsValues = sigmaSqsInput.get().getValues(); // these are all the unique sigmas
		Double[] allMusValues = meansInput.get().getValues(); // these are all the unique means

		int i = 0;
		for (Integer assignment: normalAssignments) {
			mus[i] = allMusValues[assignment];
			sigmaSqs[i] = allSigmaSqsValues[assignment];
			i++;
		}

		logP = MVNUtils.getSampleMultipleNormalLogLk(oneTraitDataArr, mus, sigmaSqs);
		
		return logP;
	}
	
	private void populateSampleData(boolean updateTipValues) {
		if (updateTipValues) {
			RealParameter oneTraitValues = oneTraitInput.get();

			/* This is the safest way to do it when there is
			a PhyloTMatrix whose species orders might not match
			the user input. But under WN there is no tree traversal
			so we can safely just grab values (see below) */
			// String[] spNames = oneTraitValues.getKeys();
			// int i = 0;
			// for (String spName: spNames) {
			//     oneTraitDataArr[i] = oneTraitValues.getValue(spName);
			// }
			oneTraitDataArr = oneTraitValues.getValues(); // this one assumes each species' values stay in the same position in the array
		}
	}
	
	@Override
	public boolean requiresRecalculation() {
		boolean dirty = true;
		return dirty;
	}
	
	@Override
	public void store() {
		System.arraycopy(oneTraitDataArr, 0, storedOneTraitDataArr, 0, oneTraitDataArr.length);
	}
	
	@Override
	public void restore() {		
		Double[] arrTmp;
		
		arrTmp = oneTraitDataArr;
		oneTraitDataArr = storedOneTraitDataArr;
		storedOneTraitDataArr = arrTmp;
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

