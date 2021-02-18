package testdrivers;

import beast.core.parameter.RealParameter;
import com.sun.deploy.util.GeneralUtil;
import contraband.math.MatrixUtilsContra;
import contraband.utils.GeneralUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import java.util.Arrays;

public class MatrixTestDriver {
    public static void main(String[] args) {
        int dim = 4;
        double[] mat1 = new double[]{1.52646482105748, -0.741938436577513, -1.2767596137151, -0.559436960260477, -1.14774568300998, -0.814191423611459, 2.74621874255422, 0.22965346979477, 2.1934018197927, 2.12624356160548, -1.28053456359236, -3.37349249083156, 0.720744524658857, 2.28811004932521, 0.697581307477588, -1.66075516057852};
        double[] mat2 = new double[]{1.65029309980434, 0.925771872907844, 0.692120613117935, 4.45175305038059, 1.8860078047097, 1.00084359136861, 0.64159191901492, 0.702641290641695, 4.17043966039917, 1.19662343659241, 1.38532573790503, 1.81955172948874, 1.58959697181671, 1.31194106031498, 1.48317172562709, 5.01717162678813};
        double[] mat3 = new double[]{1.52646482105748, -1.14774568300998, 2.1934018197927, 0.720744524658857, -0.741938436577513, -0.814191423611459, 2.12624356160548, 2.28811004932521, -1.2767596137151, 2.74621874255422, -1.28053456359236, 0.697581307477588, -0.559436960260477, 0.22965346979477, -3.37349249083156, -1.66075516057852};
        double[] mat4 = new double[]{7.61230092453671, 10.1422136645288, 0.519337234936284, 2.28298878064379, 3.06439583655447, 6.63349489194757, 0.11880383361131, 2.43755290750414, 0.699199952280018, 23.7734126473281, 8.21527396286738, 4.8234549203508, 1.67856508865952, 0.724961219821125, 1.24534633845848, 1.22914777835831};
        double[] resMat = new double[dim * dim];

        for (int i = 0; i < dim; i++) {
            for (int l = 0; l < dim; l++) {
                double sum2 = 0.0;
                for (int k = 0; k < dim; k++) {
                    double sum1 = 0;
                    for (int j = 0; j < dim; j++) {
                        sum1 += MatrixUtilsContra.getMatrixEntry(mat1, i, j, dim) * MatrixUtilsContra.getMatrixEntry(mat2, j, k, dim);
                    }
                    sum2 += sum1 * MatrixUtilsContra.getMatrixEntry(mat1, l, k, dim);
                }
                MatrixUtilsContra.setMatrixEntry(resMat, i, l, sum2 + MatrixUtilsContra.getMatrixEntry(mat4, i,l,dim), dim);
            }
        }

        System.out.println(Arrays.toString(resMat));



    }

}
