package test;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import contraband.valuewrappers.OneValueContTraits;

public class OneValueContTraitsTest {

	final static double EPSILON = 1e-4;
	private Double[] spValues, sp4, sp4Expected;
	private Double sp2, sp3, sp4OneValue;
	private Double[] trait1, trait1b, trait2, trait1Expected, trait1Expected2, trait2Expected;
	
	@Before
	public void setUp() throws Exception {
		RealParameter oneTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9, 0.0 });
		RealParameter twoTraitValues = new RealParameter(new Double[] { 4.1, 4.5, 5.9, 0.0, 5.1, 5.5, 6.9, 1.0 });
//		RealParameter twoTraitValues = new RealParameter(new Double[] { 4.1, 5.1, 4.5, 5.5, 5.9, 6.9, 0.0, 1.0 });
		String spNames = "sp1,sp2,sp3,sp4";
		String[] spNameInNewickOrder = { "sp1", "sp2", "sp3", "sp4" };
		String[] spNameInNewickOrder2 = { "sp4", "sp2", "sp3", "sp1" };
		
		OneValueContTraits oneTrait = new OneValueContTraits();
		oneTrait.initByName("nTraits", 1, "spNames", spNames, "traitValues", oneTraitValues);
		
		OneValueContTraits twoTrait = new OneValueContTraits();
		twoTrait.initByName("nTraits", 2, "spNames", spNames, "traitValues", twoTraitValues);
		
		spValues = oneTrait.getSpValues("sp3");
		System.out.println("spValues=" + Arrays.toString(spValues));
		sp3 = spValues[0];
		sp2 = oneTrait.getSpValue("sp2", 0);
		
		sp4 = twoTrait.getSpValues("sp4");
		System.out.println("sp4=" + Arrays.toString(sp4));
		sp4OneValue = twoTrait.getSpValue("sp4", 1);
	
		sp4Expected = new Double[] { 0.0, 1.0 };
		
		trait1 = twoTrait.getTraitValues(0, spNameInNewickOrder);
		trait1Expected = new Double[] { 4.1, 4.5, 5.9, 0.0 };
		
		trait1b = twoTrait.getTraitValues(0, spNameInNewickOrder2);
		trait1Expected2 = new Double[] { 0.0, 4.5, 5.9, 4.1 };
		
		trait2 = twoTrait.getTraitValues(1, spNameInNewickOrder);
		trait2Expected = new Double[] { 5.1, 5.5, 6.9, 1.0 };
	}

	@Test
	public void test1() {
		Assert.assertEquals(5.9, sp3, 0.0);
		Assert.assertEquals(4.5, sp2, 0.0);
		Assert.assertArrayEquals(sp4Expected, sp4);
	}
	
	@Test
	public void test2() {
		Assert.assertArrayEquals(trait1Expected, trait1);
		Assert.assertArrayEquals(trait1Expected2, trait1b);
		Assert.assertArrayEquals(trait2Expected, trait2);
		Assert.assertEquals(1.0, sp4OneValue, 0.0);
	}
}
