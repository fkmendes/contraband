package testdrivers;

import beast.math.matrixalgebra.IllegalDimension;
import contraband.math.LUDecompositionForArray;
import contraband.math.MatrixUtilsContra;
import contraband.utils.GeneralUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;

public class CholeskyDecompositionTestDriver {
    public static void main(String[] args) {

        double[] traitMatArr = new double[]{2.70253291718979, -0.675893270506343, -0.636978767604558, -0.424476517740056, 3.18641293667475, 1.10418672302061, -1.07756976045724, -0.375547120541696, -1.23575130814614, 5.23980774469341, -0.0470259195430953, 3.84137554540166, -4.07887622176622, 0.0263105203964659, -0.88906482129797, -1.5066823865707, -0.21433013674153, 1.80182665630146, -0.454983892686403, 2.27790911714707, 0.300388081875372, -0.125715794750515, -5.02700675875332, 2.93373818651261};
        RealMatrix traitMat = new Array2DRowRealMatrix( new double[][]{
                {2.70253291718979, -0.675893270506343, -0.636978767604558, -0.424476517740056},
                {3.18641293667475, 1.10418672302061, -1.07756976045724, -0.375547120541696},
                {-1.23575130814614, 5.23980774469341, -0.0470259195430953, 3.84137554540166},
                {-4.07887622176622, 0.0263105203964659, -0.88906482129797, -1.5066823865707},
                {-0.21433013674153, 1.80182665630146, -0.454983892686403, 2.27790911714707},
                {0.300388081875372, -0.125715794750515, -5.02700675875332, 2.93373818651261}
        });

        double[][] traitRateArr = new double[][]{
                {9.91602664865309, -1.23976549069721, -7.02552856538092, -0.430065136898062},
                {-1.23976549069721, 6.93125372300264, 4.84141921174912, -0.152690375307574},
                {-7.02552856538092, 4.84141921174912, 12.568539468804, -0.352299539678836},
                {-0.430065136898062, -0.152690375307574, -0.352299539678836, 1.63688409232851}};
        RealMatrix traitRateRM = new Array2DRowRealMatrix(traitRateArr);

        // apache commons
        CholeskyDecomposition matChol = new CholeskyDecomposition(traitRateRM);
        System.out.println("apache commons");
        RealMatrix upperMat = matChol.getLT();
        System.out.println("upper matrix");
        GeneralUtils.displayRealMatrix(upperMat);
        LUDecomposition upperMatLUD = new LUDecomposition(upperMat);
        RealMatrix dataTransformMat = upperMatLUD.getSolver().getInverse();
        RealMatrix transformedTraitsMat = traitMat.multiply(dataTransformMat);
        System.out.println("Transformed matrix");
        GeneralUtils.displayRealMatrix(transformedTraitsMat);

        // beast class
        double[][] lowerMat2DArr;
        try {
            lowerMat2DArr = (new beast.math.matrixalgebra.CholeskyDecomposition(traitRateArr)).getL();
            // caution: this returns the lower triangular form
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Numerical exception in WishartDistribution");
        }
        RealMatrix lowerRM = new Array2DRowRealMatrix(lowerMat2DArr);
        RealMatrix upperRM = lowerRM.transpose();
        System.out.println("beast class");
        GeneralUtils.displayRealMatrix(upperRM);

        double[] lowerMatArr = new double[16];
        double[] upperMatArr = new double[16];
        for(int i = 0; i < 4; i ++){
            System.arraycopy(lowerMat2DArr[i], 0, lowerMatArr, i * 4, 4);
        }
        MatrixUtilsContra.matrixTranspose(lowerMatArr, 4, upperMatArr);
        System.out.println("upper matrix 1D array" + "\n" + Arrays.toString(upperMatArr));

        double[] lu = new double [16];
        int[] pivot = new int[4];
        boolean[] evenSingular = new boolean[2];
        double[] identityMatrix = new double[16];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < 4; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, 4);
        }
        LUDecompositionForArray.ArrayLUDecomposition(upperMatArr, lu, pivot, evenSingular, 4);
        double[] dataTransformMatArr = new double[16];
        try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evenSingular[1], 4, dataTransformMatArr);
        } catch (RuntimeException e) {
            System.out.println("Singular matrix.");
        }

        double[] transformedTraitsMatArr = new double[6 * 4];
        MatrixUtilsContra.matrixMultiply(traitMatArr, dataTransformMatArr, 6, 4, transformedTraitsMatArr);
        System.out.println("Transformed 1D array" + "\n" + Arrays.toString(transformedTraitsMatArr));

    }
}
