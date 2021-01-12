package contraband.prunelikelihood;

import beast.core.*;
import beast.core.parameter.RealParameter;

import java.io.PrintStream;

public class ThresholdDifference extends CalculationNode implements Function, Loggable {
    //final public Input<RealParameter> thresholdsInput = new Input<>("threshold", "Thresholds for mapping ordered discrete traits.", Input.Validate.REQUIRED);
    final public Input<OrderedDiscreteTraits> traitsInput = new Input<>("trait","Object for ordered discrete traits", Input.Validate.REQUIRED);

    boolean needsRecompute = true;
    private double[] difference;
    private double[] storedDifference;
    private int diffDim;
    private int nrOfTraits;
    private OrderedDiscreteTraits traits;

    @Override
    public void initAndValidate() {
        traits = traitsInput.get();
        nrOfTraits = traits.getLiabilityNr();
        int nrOfThresholds = traits.getThresholds().length;
        diffDim = nrOfThresholds - nrOfTraits;
        difference = new double[diffDim];
        storedDifference = new double[diffDim];
    }

    @Override
    public int getDimension() {
        return diffDim;
    }



    @Override
    public double getArrayValue(int dim) {
        if (needsRecompute) {
            compute();
        }
        return difference[dim];
    }

    @Override
    public double[] getDoubleValues() {
        if (needsRecompute) {
            compute();
        }
        return difference;
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
        int index = 0;
        double[] thresholds = traits.getThresholds();
        for(int i = 0; i < nrOfTraits; i++){
            double[] thresholdsForTrait = traitsInput.get().getThresholdsForTrait(i, thresholds);
            for(int j = 0; j < thresholdsForTrait.length - 1; j++){
                difference[index] = thresholdsForTrait[j+1] - thresholdsForTrait[j];
                index ++;
            }
        }
        needsRecompute = false;
    }



    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
        System.arraycopy(difference, 0, storedDifference, 0, diffDim);
        super.store();
    }

    @Override
    public void restore() {
        System.arraycopy(storedDifference, 0, difference, 0, diffDim);
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }

    /**
     * Loggable interface implementation follows
     */
    @Override
    public void init(PrintStream out) {
        out.print("Diff(" + ((BEASTObject) traits.thresholdsInput.get()).getID() + ")\t");
    }

    @Override
    public void log(long sampleNr, PrintStream out) {
        //double res = variableInput.get().getValue() * variableInput.get().getValue();

        //out.print(res + "\t");
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }
}
