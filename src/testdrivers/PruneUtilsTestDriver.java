package testdrivers;

import contraband.GeneralUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class PruneUtilsTestDriver {
    public static void main(String[] args) {
        // initialize bifurcating tree
        // set branch lengths to 1.0 for simplicity

        Double[][] evolRateMat2DArray; // use a simple example for four traits
        RealMatrix evolRateMat;
        RealMatrix aMat, cMat, eMat; // initialize these with 0's
        RealVector bVec, dVec; // initialize these with 0's
        double f;

        // block for matrix A
        System.out.println("Printing matrix A:");
        // write code to initialize aMat from PruneUtils
        GeneralUtils.displayRealMatrix(aMat);

        // block for matrix E
        System.out.println("Printing matrix E:");
        // write code to initialize eMat from PruneUtils
        GeneralUtils.displayRealMatrix(eMat);

        // do the same for the other matrices
        // and the vectors

        // finally, convert this testdriver into
        // a unit test for all components (A,B,C...)
    }
}
