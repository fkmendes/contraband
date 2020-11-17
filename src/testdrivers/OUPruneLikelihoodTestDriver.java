package testdrivers;

import contraband.utils.GeneralUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class OUPruneLikelihoodTestDriver {
    public static void main(String[] args) {
        double [][] h = new double[4][4];
        h[0][0] = 10; h[0][1] = 11; h[0][2] = 12; h[0][3] = 13;
        h[0][0] = 14; h[0][1] = 15; h[0][2] = 16; h[0][3] = 17;
        h[0][0] = 18; h[0][1] = 19; h[0][2] = 20; h[0][3] = 21;
        h[0][0] = 22; h[0][1] = 23; h[0][2] = 24; h[0][3] = 25;

        RealMatrix H = new Array2DRowRealMatrix(h);

        EigenDecomposition decompH = new EigenDecomposition(H);
        RealMatrix V = decompH.getV();
        RealMatrix VT = decompH.getVT();

        GeneralUtils.displayRealMatrix(V);

        GeneralUtils.displayRealMatrix(VT);
    }
}
