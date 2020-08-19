package testdrivers;

import java.util.Arrays;

import beast.core.parameter.RealParameter;
import contraband.valuewrappers.OneValueContTraits;

public class OneValueContTraitsTestDriver {

	public static void main(String[] args) {
//		String oneTraitValues = "sp1=4.1,sp2= 4.5,sp3=5.9,sp4 =0.0 ";
//		String twoTraitValues = "sp1=4.1,sp2=4.5 ,sp4=0.0,sp3=5.9 |sp1=5.1, sp2=5.5,sp3=6.9,sp4=1.0";
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
