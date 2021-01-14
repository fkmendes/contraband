package contraband.prunelikelihood;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.math.distributions.WishartDistribution;

public class InverseMatrixGibbsOperator extends Operator {
    final public Input<Integer> dimensionInput = new Input<>("dimension", "The number of rows (columns) of the inverse matrix, representing the number of characters.", Input.Validate.REQUIRED);
    final public Input<RealParameter> valuesInput = new Input<>("values", "An array of values in the inverse matrix, representing either trait correlation or variance-covariance.", Input.Validate.REQUIRED);
    final public Input<WishartDistribution> wishartDistributionInput = new Input<>("distr", "The wishart distribution to randomly draw samples from." , Input.Validate.REQUIRED);

    private WishartDistribution wishartDistribution;
    private int dim;
    private double[][] proposal;


    @Override
    public void initAndValidate() {
        wishartDistribution = wishartDistributionInput.get();
        dim = dimensionInput.get();
        if(wishartDistribution.df.get() < dim) {
            throw new RuntimeException("InverseMatrixGibbsOperator:: Small n for large df.");
        }
        valuesInput.get().setDimension(dim * dim);
        proposal = new double[dim][dim];
    }


    @Override
    public double proposal() {
        proposal = wishartDistribution.nextWishart();
        for(int i = 0; i < dim; i ++){
            for (int j = 0; j < dim; j ++){
                valuesInput.get().setValue(i * dim + j, proposal[i][j]);
            }
        }

        return Double.POSITIVE_INFINITY;
    }

}
