package test;

import java.util.ArrayList;
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
	private double[] sp4, sp4Expected;
	private double sp3;
	private double[] trait1, trait1b, trait2, trait1Expected, trait1Expected2, trait2Expected;
	
	@Before
	public void setUp() throws Exception {
		String[] spNames = new String[] { "sp1", "sp2", "sp3", "sp4" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		String oneTraitValues = "sp1=4.1,sp2= 4.5,sp3=5.9,sp4 =0.0 ";
		String twoTraitValues = "sp1=4.1,sp2=4.5 ,sp4=0.0,sp3=5.9 |sp1=5.1, sp2=5.5,sp3=6.9,sp4=1.0";
		List<String> spNameInNewickOrder = Arrays.asList("sp1", "sp2", "sp3", "sp4");
		List<String> spNameInNewickOrder2 = Arrays.asList("sp4", "sp2", "sp3", "sp1");
		
		OneValueContTraits oneTrait = new OneValueContTraits();
		oneTrait.initByName("nTraits", 1, "taxa", taxonSet, "traitValues", oneTraitValues);
		
		OneValueContTraits twoTrait = new OneValueContTraits();
		twoTrait.initByName("nTraits", 2, "taxa", taxonSet, "traitValues", twoTraitValues);
		
		double[] spValues = oneTrait.getSpValues("sp3");
		sp3 = spValues[0];
		
		double[] sp4 = twoTrait.getSpValues("sp4");
		double[] sp4Expected = { 0.0, 1.0 };
		
		double[] trait1 = twoTrait.getTraitValues(0, spNameInNewickOrder);
		double[] trait1Expected = { 4.1, 4.5, 5.9, 0.0 };
		
		double[] trait1b = twoTrait.getTraitValues(0, spNameInNewickOrder2);
		double[] trait1Expected2 = { 0.0, 4.5, 5.9, 4.1 };
		
		double[] trait2 = twoTrait.getTraitValues(1, spNameInNewickOrder);
		double[] trait2Expected = { 5.1, 5.5, 6.9, 1.0 };
	}

	@Test
	public void test1() {
		Assert.assertEquals(5.9, sp3, 0.0);
		Assert.assertArrayEquals(sp4Expected, sp4, 0.0);
	}
	
	@Test
	public void test2() {
		Assert.assertArrayEquals(trait1Expected, trait1, 0.0);
		Assert.assertArrayEquals(trait1Expected2, trait1b, 0.0);
		Assert.assertArrayEquals(trait2Expected, trait2, 0.0);
	}
}
