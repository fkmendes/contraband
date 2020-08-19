package contraband.utils;

import contraband.math.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class GeneralUtils {
	
	public static void displayRealMatrix(RealMatrix realMat) {
		int rows = realMat.getRowDimension();
		int cols = realMat.getColumnDimension();
		
		for (int i = 0; i < rows; i++) {
		    for (int j = 0; j < cols; j++) {
		        System.out.print(realMat.getEntry(i, j) + " ");
		    }
		    
		    System.out.println();
		}
	}
	
	public static void display2DArray(double[][] mat) {
		int rows = mat.length;
		int cols = mat[0].length;
		
		for (int i = 0; i < rows; i++) {
		    for (int j = 0; j < cols; j++) {
		        System.out.print(mat[i][j] + " ");
		    }
		    
		    System.out.println();
		}
	}

	public static void displayRealVector(RealVector realVec) {
		int dim = realVec.getDimension();
		
		for (int i = 0; i < dim; i++) { 
			System.out.println(realVec.getEntry(i));
		}
	}

}
