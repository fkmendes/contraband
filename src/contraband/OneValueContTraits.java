package contraband;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class OneValueContTraits extends BEASTObject {

	final public Input<Integer> nTraitsInput = new Input<>("nTraits", "contains number of traits.", Validate.REQUIRED);
	final public Input<String> spNamesInput = new Input<>("spNames", "specifies species names separated by comma (e.g., sp1,sp2,sp3), while specifying which real parameter refers to which species in traitValues.", Validate.REQUIRED);
	final public Input<RealParameter> traitInput = new Input<>("traitValues", "quantitative trait values, one per species (if fixed, then don't operate on them).", Validate.REQUIRED);
	
	Integer nTraits, nSpp;
	String[] spNames;
	Double[] traitValues, thisSpTraitValues;
	Map<String, Double[]> spValuesMap = new HashMap<>(); 
	
	@Override
	public void initAndValidate() {
		// Getting inputs
		nTraits = nTraitsInput.get();
		spNames = spNamesInput.get().replace("\\s+", "").split(",");
		nSpp = spNames.length;
		traitValues = traitInput.get().getValues();
		
		populateSpValuesMap();
		checkAllSpHaveValues();
	}

	private void checkAllSpHaveValues() {
		;
	}
	
	private void populateSpValuesMap() {
		
		// Looping over species
		int ithSpp = 0;
		for (String spName: spNames) {
			thisSpTraitValues = Arrays.copyOfRange(traitValues, ithSpp*nTraits, ithSpp*nTraits+nTraits);
			spValuesMap.put(spName, thisSpTraitValues);
			ithSpp++;
		}
	}
	
//	private void populateSpValuesMap() {
//		String spName;
//		Double value;
//		String[] traitValueStrings = traitValueString.split("\\|");
//		
//		// Looping over traits
//		for (String oneTraitString: traitValueStrings) {
//			String[] strPairs = oneTraitString.split(","); // 0: sp name, 1: that trait value
//			
//			// Looping over species
//			for (String strPair: strPairs) {
//				String[] spAndValue = strPair.split("=");
//				spName = spAndValue[0];
//				value = Double.parseDouble(spAndValue[1]);
//				
//				// System.out.println("Parsing species " + spName + " with value " + value);
//				
//				if (!spValuesMap.containsKey(spAndValue[0])) {
//					spValuesMap.put(spName, new ArrayList<Double>());
//				}
//				
//				spValuesMap.get(spName).add(value);
//			}
//		}
//	}
	
	/*
	 * One species, one trait (need to provide trait index)
	 */
	public double getSpValue(String spName, int idxOfTraitToReturn) {
		double spValue = spValuesMap.get(spName)[idxOfTraitToReturn];

		return spValue;
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
	 * One trait, all species (order provided by spNamesInput) 
	 */
	public double[] getTraitValues(int traitIdx, String[] strings) {
//		nSpp = strings.length;
		double[] traitValues = new double[nSpp]; // used by getter (same trait, different species)
		
		int ithSpp=0;
		for (String spName: strings) {
			traitValues[ithSpp] = getSpValues(spName)[traitIdx];
			ithSpp++;
		}
		return traitValues;
	}
}
