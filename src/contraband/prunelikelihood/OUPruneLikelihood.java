package contraband.prunelikelihood;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.utils.PruneLikelihoodUtils;
import org.apache.commons.math3.linear.*;

import java.util.Arrays;

@Description("This class implements likelihood for continuous traits under Ornstein–Uhlenbeck process.\n" +
        "The calculation uses Venelin's PCM likelihood.")


public class OUPruneLikelihood extends PruneLikelihoodProcess {
    final public Input<RealParameter> sigmaValuesInput = new Input<>("sigma","Trait rate matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmaeValuesInput = new Input<>("sigmae","Population error..", Input.Validate.OPTIONAL);

    final public Input<RealParameter> alphaInput = new Input<>("alpha","An array of (nTraits * (nTraits - 1) / 2) elements, representing selection strength, off-diagonal elements in Alpha matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> thetaInput = new Input<>("theta","An array of nTraits elements, representing optimum trait values, elements in Theta vector.", Input.Validate.REQUIRED);

    private RealMatrix sigmaRM;
    private RealMatrix sigmaERM;
    private RealMatrix alphaRM;
    private RealVector thetaVec;

    private double [] sigma;
    private double [] sigmae;
    private double[] alpha;
    private double[] theta;



    @Override
    public void initAndValidate() {
        super.initAndValidate();
        int nTraits = getNTraits();

        sigmaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        populateRealMatrix(sigmaRM, sigmaValuesInput.get().getDoubleValues());

        alphaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        populateRealMatrix(alphaRM, alphaInput.get().getDoubleValues());

        thetaVec = new ArrayRealVector(new double[nTraits]);
        populateRealVector(thetaVec, thetaInput.get().getDoubleValues());

        if (sigmaeValuesInput.get() != null) {
            sigmaERM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
            populateRealMatrix(sigmaERM, sigmaeValuesInput.get().getDoubleValues());
        }

    }

    @Override
    public double calculateLogP() {

        return logP;
    }

    private void populateRealMatrix(RealMatrix rm, double[] values) {
        int k = 0;
        for (int i = 0; i < rm.getColumnDimension(); i++) {
            for (int j = i; j < rm.getColumnDimension(); j++) {
                rm.setEntry(i, j, values[k]);
                k++;
            }
        }
    }

    private void populateRealVector(RealVector vec, double[] values) {
        for (int i = 0; i < values.length; i ++) {
            vec.setEntry(i, values[i]);
        }
    }

}
