package contraband.math;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;

import beast.base.inference.distribution.ParametricDistribution;
import beast.base.core.Description;
import beast.base.core.Input;
import org.apache.commons.math3.util.FastMath;


@Description("Uniform distribution over a given interval (including lower and upper values)")
public class UniformSoftBounds extends ParametricDistribution {
    final public Input<Double> lowerInput = new Input<>("lower", "lower bound on the interval, default 0", 0.0);
    final public Input<Double> upperInput = new Input<>("upper", "lower bound on the interval, default 1", 1.0);

    UniformImpl distr = new UniformImpl();

    double _lower, _upper, density;
    private double theta1, theta2;
    private double constant1, constant2, constant3;

    private boolean infiniteSupport;

    @Override
    public void initAndValidate() {
        _lower = lowerInput.get();
        _upper = upperInput.get();
        if (_lower >= _upper) {
            throw new IllegalArgumentException("Upper value should be higher than lower value");
        }

        double diff = _upper - _lower;
        theta1 = 0.95 * _lower / (0.025 * diff);
        theta2 = 0.95 / (0.025 * diff);

        constant1 = theta1 / _lower;
        constant2 = theta1 - 1.0;
        constant3 = theta2 * _upper;



        distr.setBounds(_lower, _upper);
        infiniteSupport = Double.isInfinite(_lower) || Double.isInfinite(_upper);
        if (infiniteSupport) {
            density = 1.0;
        } else {
            density = 0.95 / (_upper - _lower);
        }
    }


    class UniformImpl implements ContinuousDistribution {
        private double lower;
        private double upper;

        public void setBounds(final double lower, final double upper) {
            this.lower = lower;
            this.upper = upper;
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
            if(x <= lower) {
                return 0.025 * constant1 * Math.pow(1 / lower, theta1 - 1) * Math.pow(x, theta1) / (theta1);
            } else if (x > upper) {
                 return -0.025 * theta2 * Math.exp(-theta2 * x + constant3);
            } else {
                return (0.95 * x - lower) / (upper - lower);
            }
        }

        @Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            x0 = Math.max(x0, lower);
            x1 = Math.min(x1, upper);
            if (x1 < lower || x1 > upper) {
                throw new RuntimeException("Value x (" + x1 + ") out of bounds (" + lower + "," + upper + ").");
            }
            return (x1 - x0) / (upper - lower);
        }

        @Override
        public double inverseCumulativeProbability(final double p) throws MathException {
            if (p < 0.0 || p > 1.0) {
                throw new RuntimeException("inverseCumulativeProbability::argument out of range [0...1]");
            }
            if( p == 0 ) {
                // works even when one bound is infinite
                return _lower;
            }
            if( p == 1 ) {
                // works even when one bound is infinite
                return _upper;
            }
            if( infiniteSupport ) {
                 throw new RuntimeException("Inverse Cumulative Probability for 0 < p < 1 and infinite support") ;
            }
            return (upper - lower) * p + lower;
        }

        @Override
        public double density(final double x) {
            if(x <= lower){
                return 0.025 * constant1 * FastMath.pow(x/lower, constant2);
            } else if (x > upper) {
                return 0.025 * theta2 * FastMath.exp(-theta2 * x + constant3);
            } else {
                return 0.95 / (upper - lower);
            }
        }

        @Override
        public double logDensity(final double x) {
            return Math.log(density(x));
        }
    } // class UniformImpl


    @Override
    public Distribution getDistribution() {
        return distr;
    }

    @Override
    public double density(final double x) {
       return density;
    }

}
