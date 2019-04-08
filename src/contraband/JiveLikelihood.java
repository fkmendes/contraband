package contraband;

import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class JiveLikelihood extends Distribution {

	final public Input<ManyValuesOneContTrait> oneTraitInput = new Input<>("sampleData", "TWO OR MORE continuous data values for ONE trait, from many species.", Validate.REQUIRED);
	final public Input<RealParameter> logSigmasqsInput = new Input<>("logSigmasqs", "log(Sigma^2), the log-variances of the normal densities, one per species.", Validate.REQUIRED);
	final public Input<RealParameter> meansInput = new Input<>("mus", "mus, the means of the normal densities, one per species.", Validate.REQUIRED);
	
	private ManyValuesOneContTrait sampleData;
	private Double[] logSigmasqs, mus; // in the order of spNames from sampleData
	private int nSpp;
	private String[] spNames;
	
	@Override
	public void initAndValidate() {
		sampleData = oneTraitInput.get();
		nSpp = sampleData.getNSpp();
		spNames = sampleData.getSpNames().toArray(new String[nSpp]);
		
		logSigmasqs = logSigmasqsInput.get().getValues();
		mus = meansInput.get().getValues();
	}
	
	@Override
	public double calculateLogP() {	
		logSigmasqs = logSigmasqsInput.get().getValues();
		mus = meansInput.get().getValues();
		
		double logLik = 0.0;
		
		int i = 0;
		for (String spName: spNames) {
			double thisSpLogLik = MVNUtils.getSampleNormalLogLk(sampleData.getSample(spName), mus[i], logSigmasqs[i]);
			logLik += thisSpLogLik;
			i++;
		}

		return logLik;
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
