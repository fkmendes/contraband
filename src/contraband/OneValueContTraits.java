package contraband;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class OneValueContTraits extends CalculationNode {

	final public Input<Integer> nTraitsInput = new Input<>("nTraits", "contains number of traits.", Validate.REQUIRED);
	final public Input<String> spNamesInput = new Input<>("spNames", "specifies species names separated by comma (e.g., sp1,sp2,sp3), while specifying which real parameter refers to which species in traitValues.", Validate.REQUIRED);
	final public Input<RealParameter> traitInput = new Input<>("traitValues", "quantitative trait values, one per species (if fixed, then don't operate on them).", Validate.REQUIRED);
	
	Integer nTraits, nSpp;
	String[] spNames;
	Double[] traitValues, thisSpTraitValues;
	Map<String, Double[]> spValuesMap;
	
	// stored stuff
	Map<String, Double[]> storedSpValuesMap;
	
	@Override
	public void initAndValidate() {
		// Getting inputs
		nTraits = nTraitsInput.get();
		thisSpTraitValues = new Double[nTraits];
		spNames = spNamesInput.get().replace("\\s+", "").split(",");
		nSpp = spNames.length;
		spValuesMap = new HashMap<String, Double[]>();
		
		populateSpValuesMap();
		checkAllSpHaveValues();
		
		storedSpValuesMap = new HashMap<String, Double[]>();
	}

	private void checkAllSpHaveValues() {
		;
	}
	
	private void populateSpValuesMap() {
		traitValues = traitInput.get().getValues();
		
		// Looping over jth traits
		for (int j=0; j < nTraits; ++j) {
		
			// Looping over ith species
			int i=0;
			for (String spName: spNames) {
				if (!spValuesMap.containsKey(spName)) {
					spValuesMap.put(spName, new Double[nTraits]);
				}
				
				spValuesMap.get(spName)[j] = traitValues[(j*nSpp) + i];
				++i;
			}
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
	
	// getters
	
	/*
	 * One species, one trait (need to provide trait index)
	 */
	public Double getSpValue(String spName, int idxOfTraitToReturn) {	
//		if (traitInput.isDirty()) {
//			populateSpValuesMap();
//		}

		Double spValue = spValuesMap.get(spName)[idxOfTraitToReturn];

		return spValue;
	}
	
	/*
	 * One species, all traits (same order as in string)
	 */
	public Double[] getSpValues(String spName) {	
//		if (traitInput.isDirty()) {
//			populateSpValuesMap();
//		}

		return spValuesMap.get(spName);
	}
	
	/*
	 * One trait, all species (order provided by spNamesInput) 
	 */
	public Double[] getTraitValues(int traitIdx, String[] strings) {
//		if (traitInput.isDirty()) {
//			populateSpValuesMap();
//		}

		Double[] traitValues = new Double[nSpp]; // used by getter (same trait, different species)
		
		int ithSpp=0;
		for (String spName: strings) {
			traitValues[ithSpp] = getSpValues(spName)[traitIdx];
			ithSpp++;
		}
		
		System.out.println("Species names for trait values below:" + Arrays.toString(strings));
		System.out.println("traitValues inside OneValueContTraits=" + Arrays.toString(traitValues));
		
		return traitValues;
	}
		
	public Integer getNTraits() {
		return nTraits;
	}
	
	// caching
	@Override
	protected boolean requiresRecalculation() {
		return true;
	}
	
	@Override
	public void store() {
		for (String spName: spValuesMap.keySet()) {
			storedSpValuesMap.put(spName, spValuesMap.get(spName).clone());
		}
		
		super.store();
	}
	
	@Override
	public void restore() {
		Map<String, Double[]> mapTmp;
		
		mapTmp = spValuesMap;
		spValuesMap = storedSpValuesMap;
		storedSpValuesMap = mapTmp;
		
		super.restore();
	}
}
