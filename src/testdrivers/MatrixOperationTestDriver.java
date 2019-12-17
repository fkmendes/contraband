package testdrivers;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import contraband.MatrixOperation;
import contraband.GeneralUtils;


public class MatrixOperationTestDriver {

	public static void main(String[] args) {

		RealMatrix aMat = new Array2DRowRealMatrix(new double [][]
				{{1.0, 2.0},{3.0, 4.0}}
		);

		RealMatrix bMat = new Array2DRowRealMatrix(new double [][]
				{{5.0, 6.0},{7.0, 8.0}}
		);
		RealMatrix resRM = new Array2DRowRealMatrix(new double [2][2]);
		// block for matrix add
		System.out.println("Printing matrix aMat:");
		GeneralUtils.displayRealMatrix(aMat);
		System.out.println("Printing matrix bMat:");
		GeneralUtils.displayRealMatrix(bMat);
		MatrixOperation.matrixAdd(aMat, bMat, resRM);
		System.out.println("Printing matrix aMat + bMat:");
		GeneralUtils.displayRealMatrix(resRM);


		// block for matrix subtract
		System.out.println("Printing matrix A before subtract:");
		GeneralUtils.displayRealMatrix(aMat);
		MatrixOperation.matrixSubtract(aMat, bMat);
		System.out.println("Printing matrix A after aMat - bMat:");
		GeneralUtils.displayRealMatrix(aMat);

		// block for matrix scalar multiply
		System.out.println("Printing matrix A before multiply:");
		GeneralUtils.displayRealMatrix(aMat);
		MatrixOperation.matrixScalarMultiply(aMat, 2);
		System.out.println("Printing matrix A after aMat * 2:");
		GeneralUtils.displayRealMatrix(aMat);

		// block for matrix scalar add
		System.out.println("Printing matrix A before scalar add:");
		GeneralUtils.displayRealMatrix(aMat);
		MatrixOperation.matrixScalarAdd(aMat, bMat, 2);
		System.out.println("Printing matrix A after aMat + 2 * bMat:");
		GeneralUtils.displayRealMatrix(aMat);

		// block for matrix transpose
		System.out.println("Printing matrix A before transpose:");
		GeneralUtils.displayRealMatrix(aMat);
		MatrixOperation.matrixTranspose(aMat);
		System.out.println("Printing matrix A after aMat^T:");
		GeneralUtils.displayRealMatrix(aMat);

	}
}
