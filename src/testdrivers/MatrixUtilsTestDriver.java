package testdrivers;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import contraband.MatrixUtils;
import contraband.GeneralUtils;


public class MatrixUtilsTestDriver {

	public static void main(String[] args) {

		RealMatrix inRM = new Array2DRowRealMatrix(new double [][]
				{{1.0, 2.0}, {3.0, 4.0}}
		);

		RealMatrix inRM2 = new Array2DRowRealMatrix(new double [][]
				{{10.0, 20.0}, {30.0, 40.0}}
		);

		RealMatrix rmToAdd = new Array2DRowRealMatrix(new double [][]
				{{5.0, 6.0}, {7.0, 8.0}}
		);

		RealMatrix resRM = new Array2DRowRealMatrix(new double [2][2]);

		RealMatrix resRM2 = new Array2DRowRealMatrix(new double [2][2]);

		// printing inputs and result for sum
		System.out.println("Printing matrix inRM:");
		GeneralUtils.displayRealMatrix(inRM);
		System.out.println("Printing matrix rmToAdd:");
		GeneralUtils.displayRealMatrix(rmToAdd);
		resRM = MatrixUtils.matrixAdd(inRM, rmToAdd, resRM);
		System.out.println("Printing matrix inRM + rmToAdd:");
		GeneralUtils.displayRealMatrix(resRM);

		// doing another sum, printing result to see if it changed
		// (checking side-effect)
		resRM2 = MatrixUtils.matrixAdd(inRM2, rmToAdd, resRM);
		System.out.println("Printing original result matrix to see if there was a side-effect:");
		GeneralUtils.displayRealMatrix(resRM);
		System.out.println("Printing second result matrix:");
		GeneralUtils.displayRealMatrix(resRM2);

		// printing result for subtraction
		resRM = MatrixUtils.matrixSubtract(inRM, rmToAdd, resRM);
		System.out.println("Printing matrix inRM - rmToAdd:");
		GeneralUtils.displayRealMatrix(resRM);

		// printing result for scalar multiply
		resRM = MatrixUtils.matrixScalarMultiply(inRM, 2, resRM);
		System.out.println("Printing matrix inRM * 2:");
		GeneralUtils.displayRealMatrix(resRM);

		// printing result for multiply
		resRM = MatrixUtils.matrixMultiply(inRM, rmToAdd, resRM);
		System.out.println("Printing matrix inRM * rmToAdd:");
		GeneralUtils.displayRealMatrix(resRM);

		// printing result for scalar add
		MatrixUtils.matrixScalarAdd(inRM, 2, resRM);
		System.out.println("Printing matrix inRM + 2:");
		GeneralUtils.displayRealMatrix(resRM);

		// printing result for transpose
		MatrixUtils.matrixTranspose(inRM, resRM);
		System.out.println("Printing matrix inRM^T:");
		GeneralUtils.displayRealMatrix(resRM);
	}
}
