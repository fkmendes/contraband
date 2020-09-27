package contraband.otherlikelihood;

import java.util.List;
import java.util.Random;

import beast.core.*;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import contraband.math.MVNUtils;
import contraband.valuewrappers.ManyValuesOneContTrait;

/**
 * @author Fabio K. Mendes
 */

@Description("JIVE model (Joint inter- and Intraspecific Variance Evolution) implementation," +
        "for analyses of continuous trait evolution in multiple individuals" +
        "and multiple species.")
@Citation(value = "Gaboriau, T. et al. (2020). A multi-platform package for" +
        "the analysis of intra- and interspecific trait evolution. Method. Ecol. Evol., 1-9",
        DOI = "10.1111/2041-210X.13458",
        year = 2020,
        firstAuthorSurname = "Gaboriau")
public class JiveLikelihoodOneTrait extends Distribution {

	final public Input<ManyValuesOneContTrait> oneTraitInput = new Input<>("sampleData", "TWO OR MORE continuous data values for ONE trait, from many species.", Validate.REQUIRED);
	final public Input<RealParameter> logSigmaSqsInput = new Input<>("logSigmaSqs", "log-sigma^2s, the log-variances of the normal densities, one per species.", Validate.REQUIRED);
	final public Input<RealParameter> meansInput = new Input<>("mus", "mus, the means of the normal densities, one per species.", Validate.REQUIRED);
	
	private ManyValuesOneContTrait sampleData;
	private int nSpp;
	private String[] spNames;
	
	@Override
	public void initAndValidate() {
		sampleData = oneTraitInput.get();
		nSpp = sampleData.getNSpp();
		spNames = sampleData.getSpNames().toArray(new String[nSpp]);
	}
	
	@Override
	public double calculateLogP() {	
		Double[] logSigmaSqs = logSigmaSqsInput.get().getValues();
		Double[] mus = meansInput.get().getValues();
		
		logP = 0.0;
		
		int i = 0;
		for (String spName: spNames) {
			double thisSpLogLik = MVNUtils.getSampleNormalLogLk(sampleData.getSample(spName), mus[i], logSigmaSqs[i]);
			logP += thisSpLogLik;
			i++;
		}

		return logP;
	}
	
	@Override
	protected boolean requiresRecalculation() {
		// TODO Auto-generated method stub
		return super.requiresRecalculation();
	}

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	
}

