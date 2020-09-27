package testdrivers;

import beast.core.parameter.RealParameter;
import contraband.otherlikelihood.JiveLikelihoodOneTrait;
import contraband.valuewrappers.ManyValuesOneContTrait;

/**
 * @author Fabio K. Mendes
 */

public class JiveLikelihoodTestDriver {

	public static void main(String[] args) {
		String oneTraitValues = "sp1=-0.65187041,2.04994431,0.70321432,0.61298809,0.04211039|sp2=0.9181909,1.7149411,-0.2674609,0.8641253,1.4442747";
		ManyValuesOneContTrait oneTrait = new ManyValuesOneContTrait();
		oneTrait.initByName("traitValues", oneTraitValues);
		
		// log-sigmasq
		Double[] logSigmasqsInput = new Double[] { 0.1, 0.15 };
		RealParameter logSigmasqs = new RealParameter(logSigmasqsInput);
				
		// mean vector
		Double[] meansInput = new Double[] { 0.1, 1.1 };
		RealParameter mus = new RealParameter(meansInput);
		
		JiveLikelihoodOneTrait jive = new JiveLikelihoodOneTrait();
		jive.initByName("sampleData", oneTrait, "logSigmaSqs", logSigmasqs, "mus", mus);
		
		/*
		 *  In R:
		 *  data1 <- c(-0.65187041,2.04994431,0.70321432,0.61298809,0.04211039)
		 *  data2 <- c(0.9181909,1.7149411,-0.2674609,0.8641253,1.4442747)
		 *  sum(dnorm(data1, 0.1, sqrt(exp(0.1)), log=TRUE), dnorm(data2, 1.1, sqrt(exp(0.15)), log=TRUE))
		 */
		System.out.println(jive.calculateLogP()); // from R: -13.13222
	}
}
