package testdrivers;

import beast.core.parameter.RealParameter;
import beast.math.distributions.WishartDistribution;
import contraband.prunelikelihood.InverseMatrixGibbsOperator;

import java.util.Arrays;

public class InverseMatrixGibbsOperatorTestDriver {
    public static void main(String[] args) {
        WishartDistribution wishartDistribution = new WishartDistribution();
        RealParameter scaleMatrix  = new RealParameter(new Double[] {
                1.0, -0.364849804288319, -0.22682370058509, 0.0132573574252411, 0.22906744642656, 0.545701377465508,
                -0.364849804288319, 1.0, -0.120181641613719, -0.0537910912563028, -0.222377622324743, -0.645977169868794,
                -0.22682370058509, -0.120181641613719, 1.0, 0.553396782108105, 0.0767992903301167, -0.151551897322574,
                0.0132573574252411, -0.0537910912563028, 0.553396782108105, 1.0, -0.0358374883713417, -0.0844533785837393,
                0.22906744642656, -0.222377622324743, 0.0767992903301167, -0.0358374883713417, 1.0, -0.248251702214294,
                0.545701377465508, -0.645977169868794, -0.151551897322574, -0.0844533785837393, -0.248251702214294, 1.0
        });
        wishartDistribution.initByName("df", 6.0, "scaleMatrix", scaleMatrix);

        InverseMatrixGibbsOperator operator = new InverseMatrixGibbsOperator();
        RealParameter inverseMatrix  = new RealParameter(new Double[] {1.0});
        operator.initByName("dimension", 6, "values", inverseMatrix, "distr", wishartDistribution, "weight", 1.0);

        operator.proposal();
        System.out.println(Arrays.toString(inverseMatrix.getDoubleValues()));
    }
}
