package test;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import contraband.JiveLikelihood;
import contraband.ManyValuesOneContTrait;

public class JiveLikelihoodTest {

	double logLik;
	
	@Before
	public void setUp() throws Exception {
		String oneTraitValues = "sp1=-0.65187041,2.04994431,0.70321432,0.61298809,0.04211039|sp2=0.9181909,1.7149411,-0.2674609,0.8641253,1.4442747";
		ManyValuesOneContTrait oneTrait = new ManyValuesOneContTrait();
		oneTrait.initByName("traitValues", oneTraitValues);
		
		// log-sigmasq
		Double[] logSigmasqsInput = new Double[] { 0.1, 0.15 };
		RealParameter logSigmasqs = new RealParameter(logSigmasqsInput);
				
		// mean vector
		Double[] meansInput = new Double[] { 0.1, 1.1 };
		RealParameter mus = new RealParameter(meansInput);
		
		JiveLikelihood jive = new JiveLikelihood();
		jive.initByName("sampleData", oneTrait, "logSigmasqs", logSigmasqs, "mus", mus);
		
		logLik = jive.calculateLogP();
	}

	@Test
	public void test() {
		assertEquals(-13.13222, logLik, 1E-5);
	}

}