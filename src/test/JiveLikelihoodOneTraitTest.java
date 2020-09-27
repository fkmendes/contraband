package test;

import static org.junit.Assert.*;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import contraband.otherlikelihood.JiveLikelihoodOneTrait;
import contraband.valuewrappers.ManyValuesOneContTrait;

/**
 * @author Fabio K. Mendes
 */

public class JiveLikelihoodOneTraitTest {

	final static double EPSILON = 1e-6;
	
	@Test
	public void testJiveLkOneTrait() {
		String oneTraitValues = "sp1=-0.65187041,2.04994431,0.70321432,0.61298809,0.04211039|sp2=0.9181909,1.7149411,-0.2674609,0.8641253,1.4442747";
		ManyValuesOneContTrait oneTrait = new ManyValuesOneContTrait();
		oneTrait.initByName("traitValues", oneTraitValues);
		
		// sigmasq
		Double[] logSigmaSqsInput = new Double[] { 0.1, 0.15 };
		RealParameter logSigmaSqs = new RealParameter(logSigmaSqsInput);
				
		// mean vector
		Double[] meansInput = new Double[] { 0.1, 1.1 };
		RealParameter mus = new RealParameter(meansInput);
		
		JiveLikelihoodOneTrait jive = new JiveLikelihoodOneTrait();
		jive.initByName("sampleData", oneTrait, "logSigmaSqs", logSigmaSqs, "mus", mus);
		double logLik = jive.calculateLogP();

		assertEquals(-13.1322153, logLik, EPSILON);
	}
}
