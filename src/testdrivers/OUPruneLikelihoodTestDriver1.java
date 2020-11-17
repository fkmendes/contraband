package testdrivers;

import contraband.utils.GeneralUtils;
import org.apache.commons.math3.linear.*;
import java.util.Arrays;

public class OUPruneLikelihoodTestDriver1 {
    public static void main(String[] args) {
        double [][] h = new double[4][4];
        h[0][0] = 10; h[0][1] = 11; h[0][2] = 12; h[0][3] = 13;
        h[1][0] = 14; h[1][1] = 15; h[1][2] = 16; h[1][3] = 17;
        h[2][0] = 18; h[2][1] = 19; h[2][2] = 20; h[2][3] = 21;
        h[3][0] = 22; h[3][1] = 23; h[3][2] = 24; h[3][3] = 25;

        RealMatrix H = new Array2DRowRealMatrix(h);
        System.out.println("Display H: ");
        GeneralUtils.displayRealMatrix(H);

        //H = V × D × VT
        EigenDecomposition decompH = new EigenDecomposition(H);

        // (1)
        System.out.println("Display V: ");
        RealMatrix V = decompH.getV();
        GeneralUtils.displayRealMatrix(V);

        // (2)
        System.out.println("Display VT: ");
        RealMatrix VT = decompH.getVT();
        GeneralUtils.displayRealMatrix(VT);

        // (3)
        System.out.println("Display D: ");
        RealMatrix D = decompH.getD();
        GeneralUtils.displayRealMatrix(D);

        // (4)
        RealVector vec1 = decompH.getEigenvector(0);
        RealVector vec2 = decompH.getEigenvector(1);
        RealVector vec3 = decompH.getEigenvector(2);
        RealVector vec4 = decompH.getEigenvector(3);
        System.out.println("Print vec1:");
        GeneralUtils.displayRealVector(vec1);
        System.out.println("Print vec2:");
        GeneralUtils.displayRealVector(vec2);
        System.out.println("Print vec3:");
        GeneralUtils.displayRealVector(vec3);
        System.out.println("Print vec4:");
        GeneralUtils.displayRealVector(vec4);

        // (5)
        double[] values = decompH.getRealEigenvalues();
        System.out.println("Display eigen values: " + "\n" + Arrays.toString(values));

        // (6)
        System.out.println("Print inverse matrix:");
        RealMatrix inverse = new LUDecomposition(V).getSolver().getInverse();
        GeneralUtils.displayRealMatrix(inverse);
    }
}
