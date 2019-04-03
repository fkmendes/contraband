package contraband;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Input.Validate;

public class ManyValuesOneContTrait extends BEASTObject {

	final public Input<String> traitInput = new Input<>("traitValues", "contains String array, each string containing 2+ values for one of the traits for all species in sp1=value1,value2|sp2=value1,value2,etc. format.", Validate.REQUIRED);

	String traitValueString;
	String[] thisSpeciesValues;
	double[] spValues; // return
	int nSpp;
	List<String> spNames;
	Map<String, List<Double>> spValuesMap = new HashMap<>(); 
	
	@Override
	public void initAndValidate() {
		// Getting inputs
		traitValueString = traitInput.get().replaceAll("\\s+","");
		spNames = new ArrayList<String>();
		
		populateSpValuesMap();
		checkAllSpHave2PlusValues();
	}

	private void populateSpValuesMap() {
		String[] eachSpeciesStuff = traitValueString.split("\\|");
		String spName;
		
		// Looping over species
		for (String oneSpeciesSamples: eachSpeciesStuff) {
			String[] strTokens = oneSpeciesSamples.split("="); // 0: sp name, 1: string with samples
			spName = strTokens[0];
			spNames.add(spName);
			thisSpeciesValues = strTokens[1].split(",");
			
			// Looping over samples within a species
			for (String oneValue: thisSpeciesValues) {
				if (!spValuesMap.containsKey(spName)) {
					spValuesMap.put(spName, new ArrayList<Double>());
				}
				
				spValuesMap.get(spName).add(Double.valueOf(oneValue));
			}
		}
	}

	/*
	 * Get all values from the one trait from a species
	 */
	public double[] getSample(String spName) {
		spValues = new double[spValuesMap.get(spName).size()]; // used by getter (different traits, same species)

		int i=0;
		for (Double spValue: spValuesMap.get(spName)) {
			spValues[i] = spValue.doubleValue();
			i++;
		}

		return spValues;
	}
	
	public int getNSpp() {
		return spValuesMap.size();
	}
	
	public List<String> getSpNames() {
		return spNames;
	}
	
	private void checkAllSpHave2PlusValues() {
		// TODO Auto-generated method stub
		;
	}
	
	
}
