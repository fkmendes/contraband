package testdrivers;

import java.util.Arrays;
import java.util.List;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import contraband.OneValueContTraits;

public class OneValueContTraitsTestDriver {

	public static void main(String[] args) {
		String oneTraitValues = "sp1=4.1,sp2= 4.5,sp3=5.9,sp4 =0.0 ";
		String twoTraitValues = "sp1=4.1,sp2=4.5 ,sp4=0.0,sp3=5.9 |sp1=5.1, sp2=5.5,sp3=6.9,sp4=1.0";
		
		OneValueContTraits oneTrait = new OneValueContTraits();
		oneTrait.initByName("nTraits", 1, "traitValues", oneTraitValues);
		
		OneValueContTraits twoTrait = new OneValueContTraits();
		twoTrait.initByName("nTraits", 2, "traitValues", twoTraitValues);
		
		double[] spValues = oneTrait.getSpValues("sp3");
		double[] spValues2 = twoTrait.getSpValues("sp4");
		
		System.out.println(Arrays.toString(spValues));
		System.out.println(Arrays.toString(spValues2));
	}

}
