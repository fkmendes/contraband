package contraband.math;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;


@Description("calculates sum of a valuable")
public class Square extends CalculationNode implements Function, Loggable {
    final public Input<RealParameter> variableInput = new Input<>("arg", "argument to be summed", Validate.REQUIRED);


    boolean needsRecompute = true;
    double square = 0;
    double storedSquare = 0;


    @Override
    public void initAndValidate() {
        RealParameter valuable = variableInput.get();
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return square;
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {

        square = variableInput.get().getValue() * variableInput.get().getValue();

        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (dim == 0) {
            return getArrayValue();
        }
        return Double.NaN;
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
        storedSquare = square;
        super.store();
    }

    @Override
    public void restore() {
        square = storedSquare;
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
        out.print("square(" + ((BEASTObject) variableInput.get()).getID() + ")\t");
    }

    @Override
    public void log(long sampleNr, PrintStream out) {
	    double res = variableInput.get().getValue() * variableInput.get().getValue();

        out.print(res + "\t");
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

} // class Sum
