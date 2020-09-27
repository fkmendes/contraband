package test;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;
import contraband.coalescent.CoalUtils;

/**
 * @author Fabio K. Mendes
 */

public class CoalUtilsTest {

	double g10to2tau3Tavare, g10to2t3000ne1000Heled;
	
	@Before
	public void setUp() throws Exception {
		g10to2tau3Tavare = CoalUtils.getTavareGij(10, 2, 3);
		g10to2t3000ne1000Heled = CoalUtils.getHeledGij(10, 2, 3000, 1000);
	}

	@Test
	public void test() {
		assertEquals(g10to2tau3Tavare, g10to2t3000ne1000Heled, 1E-8);
		assertEquals(0.12119509364945, g10to2tau3Tavare, 1E-8);
		assertEquals(0.12119509364945, g10to2t3000ne1000Heled, 1E-8);
	}
}
