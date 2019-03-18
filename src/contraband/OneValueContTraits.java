package contraband;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.TaxonSet;

public class OneValueContTraits extends BEASTObject {

	final public Input<Integer> nTraitsInput = new Input<>("nTraits", "contains number of traits.", Validate.REQUIRED);
	final public Input<String> traitInput = new Input<>("traitValues", "contains String array, each string containing values for one of the traits for all species in sp_name=value form.", Validate.REQUIRED);
	final public Input<TaxonSet> taxaInput = new Input<>("taxa", "contains list of taxa to map traits to.", Validate.REQUIRED);
	
	Integer nTraits, nSpp;
	String traitValueString;
	String[] traitValueStrings;
	Set<String> spNames;
	Map<String, List<Double>> spValuesMap = new HashMap<>(); 
	
	@Override
	public void initAndValidate() {
		// Getting inputs
		nTraits = nTraitsInput.get();
		nSpp = taxaInput.get().getTaxonCount();
		spNames = taxaInput.get().getTaxaNames();
		System.out.println(spNames);
		traitValueString = traitInput.get().replaceAll("\\s+","");
		
		populateSpValuesMap();
		checkAllSpHaveValues();
	}

	private void checkAllSpHaveValues() {
		;
	}
	
	private void populateSpValuesMap() {
		String spName;
		Double value;
		String[] traitValueStrings = traitValueString.split("\\|");
		
		// Looping over traits
		for (String oneTraitString: traitValueStrings) {
			String[] strPairs = oneTraitString.split(","); // 0: sp name, 1: that trait value
			
			// Looping over species
			for (String strPair: strPairs) {
				String[] spAndValue = strPair.split("=");
				spName = spAndValue[0];
				value = Double.parseDouble(spAndValue[1]);
				
				if (!spValuesMap.containsKey(spAndValue[0])) {
					spValuesMap.put(spName, new ArrayList<Double>());
				}
				
				spValuesMap.get(spName).add(value);
				// System.out.println(spValuesMap.get(spName));
			}
		}
	}
	
	/*
	 * One species, all traits (same order as in string)
	 */
	public double[] getSpValues(String spName) {		
		double[] spValues = new double[nTraits]; // used by getter (different traits, same species)

		
		int i=0;
		for (Double spValue: spValuesMap.get(spName)) {
			spValues[i] = spValue.doubleValue();
			i++;
		}

		return spValues;
	}
	
	/*
	 * One trait, all species (in TaxonSet order)
	 */
	public double[] getTraitValues(int traitIdx, List<String> spNamesInRightOrder) {
		double[] traitValues = new double[nSpp]; // used by getter (same trait, different species)
		
		int i=0;
		for (String spName: spNamesInRightOrder) {
			traitValues[i] = getSpValues(spName)[traitIdx];
			i++;
		}
		return traitValues;
	}
}
