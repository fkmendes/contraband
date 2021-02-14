package testdrivers;

import contraband.math.MatrixUtilsContra;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import java.util.Arrays;

public class OUModelParamsTestDriver3 {

    public static void main(String[] args) {
        int nTraits = 4;
        double branchLength = 3.0058179;


        double[] alphaMat = new double[] {-0.33856519557947, -0.077219810489012, -0.336935347137926, 0.372606734891398, -0.077219810489012, -0.189182637261135, 0.277006730470855, -0.030511296063633, -0.336935347137926, 0.277006730470855, -0.338861381957585, 0.101393849483658, 0.372606734891398, -0.030511296063633, 0.101393849483658, 0.404426904243343};
        double[][] res = new double[nTraits][nTraits];
        for (int i = 0; i < nTraits; i++) {
            for (int j = 0; j < nTraits; j++) {
                res[i][j] = MatrixUtilsContra.getMatrixEntry(alphaMat, i, j, nTraits) * (-branchLength);
            }
        }

        double[][] phiMatrix = MatrixFunctions.expm(new DoubleMatrix(res)).toArray2();

        double[] phiMat = new double[] {5.79139780244504, -0.655933465833111, 4.23336767300133, -2.23847639400205, -0.655933465833112, 2.52293429844627, -2.23165633273757, 0.372164773841502, 4.23336767300133, -2.23165633273757, 5.6575229085707, -1.6836916603114, -2.23847639400205, 0.372164773841502, -1.6836916603114, 1.08811139295935};

        double[] thetaVec = new double[] {1.82708241754003, 0.0028042328310563, 1.20255228727526, -2.14754518996932};
        double[] omega = new double[nTraits];
        double[] negativeEPhi = new double[nTraits * nTraits];
        for (int i = 0; i < nTraits; i++){
            for (int j = 0; j < nTraits; j++){
                if(i == j){
                    negativeEPhi[i * nTraits + j] = 1.0 - MatrixUtilsContra.getMatrixEntry(phiMat, i, j, nTraits);
                } else {
                    negativeEPhi[i * nTraits + j] = - MatrixUtilsContra.getMatrixEntry(phiMat, j, i, nTraits);
                }
            }
        }
        MatrixUtilsContra.matrixPreMultiply(thetaVec, negativeEPhi, nTraits, nTraits, omega);
        System.out.println(Arrays.toString(omega));
    }
}
