package test;

import java.util.Arrays;
import java.util.List;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import contraband.OneValueContTraits;

public class OneValueContTraitsTest {

	final static double EPSILON = 1e-4;
	private double[] spValues, sp4, sp4Expected;
	private double sp2, sp3, sp4OneValue;
	private double[] trait1, trait1b, trait2, trait1Expected, trait1Expected2, trait2Expected;
	
	@Before
	public void setUp() throws Exception {
		String oneTraitValues = "sp1=4.1,sp2= 4.5,sp3=5.9,sp4 =0.0 ";
		String twoTraitValues = "sp1=4.1,sp2=4.5 ,sp4=0.0,sp3=5.9 |sp1=5.1, sp2=5.5,sp3=6.9,sp4=1.0";
		String[] spNameInNewickOrder = { "sp1", "sp2", "sp3", "sp4" };
		String[] spNameInNewickOrder2 = { "sp4", "sp2", "sp3", "sp1" };
		
		OneValueContTraits oneTrait = new OneValueContTraits();
		oneTrait.initByName("nTraits", 1, "traitValues", oneTraitValues);
		
		OneValueContTraits twoTrait = new OneValueContTraits();
		twoTrait.initByName("nTraits", 2, "traitValues", twoTraitValues);
		
		spValues = oneTrait.getSpValues("sp3");
		sp3 = spValues[0];
		sp2 = oneTrait.getSpValue("sp2", 0);
		
		sp4 = twoTrait.getSpValues("sp4");
		sp4OneValue = twoTrait.getSpValue("sp4", 1);
		
		sp4Expected = new double[] { 0.0, 1.0 };
		
		trait1 = twoTrait.getTraitValues(0, spNameInNewickOrder);
		trait1Expected = new double[] { 4.1, 4.5, 5.9, 0.0 };
		
		trait1b = twoTrait.getTraitValues(0, spNameInNewickOrder2);
		trait1Expected2 = new double[] { 0.0, 4.5, 5.9, 4.1 };
		
		trait2 = twoTrait.getTraitValues(1, spNameInNewickOrder);
		trait2Expected = new double[] { 5.1, 5.5, 6.9, 1.0 };
	}

	@Test
	public void test1() {
		Assert.assertEquals(5.9, sp3, 0.0);
		Assert.assertEquals(4.5, sp2, 0.0);
		Assert.assertArrayEquals(sp4Expected, sp4, 0.0);
	}
	
	@Test
	public void test2() {
		Assert.assertArrayEquals(trait1Expected, trait1, 0.0);
		Assert.assertArrayEquals(trait1Expected2, trait1b, 0.0);
		Assert.assertArrayEquals(trait2Expected, trait2, 0.0);
		Assert.assertEquals(1.0, sp4OneValue, 0.0);
	}
}
