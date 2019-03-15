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
	private double[] sp4, sp4Expected;
	private double sp3;
	
	@Before
	public void setUp() throws Exception {
		String[] spNames = new String[] { "sp1", "sp2", "sp3", "sp4" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		String oneTraitValues = "sp1=4.1,sp2= 4.5,sp3=5.9,sp4 =0.0 ";
		String twoTraitValues = "sp1=4.1,sp2=4.5 ,sp4=0.0,sp3=5.9 |sp1=5.1, sp2=5.5,sp3=6.9,sp4=1.0";
		
		OneValueContTraits oneTrait = new OneValueContTraits();
		oneTrait.initByName("nTraits", 1, "taxa", taxonSet, "traitValues", oneTraitValues);
		
		OneValueContTraits twoTrait = new OneValueContTraits();
		twoTrait.initByName("nTraits", 2, "taxa", taxonSet, "traitValues", twoTraitValues);
		
		double[] spValues = oneTrait.getSpValues("sp3");
		sp3 = spValues[0];
		
		double[] sp4 = twoTrait.getSpValues("sp4");
		double[] sp4Expected = { 0.0, 1.0 };
	}

	@Test
	public void test() {
		Assert.assertEquals(5.9, sp3, 0.0);
		Assert.assertArrayEquals(sp4Expected, sp4, 0.0);
	}
}
