package testdrivers;

import java.util.Arrays;

import beast.base.inference.parameter.RealParameter;
import contraband.valuewrappers.OneValueContTraits;

/**
 * @author Fabio K. Mendes
 */

public class OneValueContTraitsTestDriver {

	public static void main(String[] args) {
		String spNames = "sp1,sp2,sp3,sp4";
		RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9, 0.0 });
		RealParameter twoTraitValues = new RealParameter(new Double[] { 4.1, 5.1, 4.5, 5.5, 5.9, 6.9, 0.0, 1.0 });
		
		OneValueContTraits oneTrait = new OneValueContTraits();
		oneTrait.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		OneValueContTraits twoTrait = new OneValueContTraits();
		twoTrait.initByName("nTraits", 2, "spNames", spNames, "traitValues", twoTraitValues);
		
		Double[] spValues = oneTrait.getSpValues("sp3");
		Double[] spValues2 = twoTrait.getSpValues("sp4");
		double spValue3 = twoTrait.getSpValue("sp2", 1);
		
		System.out.println(Arrays.toString(spValues));
		System.out.println(Arrays.toString(spValues2));
		System.out.println(spValue3);
	}
}
