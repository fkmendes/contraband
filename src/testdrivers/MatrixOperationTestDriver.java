package testdrivers;

import contraband.utils.GeneralUtils;
import org.apache.commons.math3.linear.*;

public class MatrixOperationTestDriver {
    public static void main(String[] args) {
        System.out.println("Printing a = ");
        RealMatrix aMat = new Array2DRowRealMatrix(new double [][]
                {{1.0, 5.0, 9.0,  13.0},
                 {2.0, 6.0, 10.0, 14.0},
                 {3.0, 7.0, 11.0, 15.0},
                 {4.0, 8.0, 12.0, 16.0}});
        GeneralUtils.displayRealMatrix(aMat);

        System.out.println("Printing b = ");
        RealMatrix bMat = new Array2DRowRealMatrix(new double [][]
                       {{4.0, 8.0, 12.0, 16.0},
                        {5.0, 7.0, 11.0, 15.0},
                        {2.0, 6.0, 10.0, 14.0},
                        {1.0, 5.0, 9.0, 13.0}});
        GeneralUtils.displayRealMatrix(bMat);

        System.out.println("Printing a + b = ");
        RealMatrix aPlusBMat = aMat.add(bMat);
        GeneralUtils.displayRealMatrix(aPlusBMat);
        System.out.println("if a == aPlusBMat " + aPlusBMat.equals(aMat));

        System.out.println("Printing a - b = ");
        RealMatrix aSubtractBMat = aMat.subtract(bMat);
        GeneralUtils.displayRealMatrix(aSubtractBMat);
        System.out.println("if a == aSubtractBMat " + aSubtractBMat.equals(aMat));

        System.out.println("Printing 0.1 * a = ");
        RealMatrix aScalerMultiply = aMat.scalarMultiply(0.1);
        GeneralUtils.displayRealMatrix(aScalerMultiply);
        System.out.println("if a == aScalerMultiply " + aScalerMultiply.equals(aMat));

        System.out.println("Printing 0.1 + a = ");
        RealMatrix aScalerAdd = aMat.scalarAdd(2.0);
        GeneralUtils.displayRealMatrix(aScalerAdd);
        System.out.println("if a == aScalerAdd " + aScalerAdd.equals(aMat));

        System.out.println("Printing a.transpose = ");
        RealMatrix aTranspose = aMat.transpose();
        GeneralUtils.displayRealMatrix(aTranspose);
        System.out.println("if a == aTranspose " + aTranspose.equals(aMat));
    }
}
