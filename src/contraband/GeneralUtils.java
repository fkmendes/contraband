package contraband;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class GeneralUtils {
	
	public static void displayRealMatrix(RealMatrix Mat) {
		
		int rows = Mat.getRowDimension();
		int cols = Mat.getColumnDimension();
		
		for (int i = 0; i < rows; i++) {
			
		    for (int j = 0; j < cols; j++) {
		    	
		        System.out.print(Mat.getEntry(i, j) + " ");
		        
		    }
		    
		    System.out.println();
		}
		
	}
	
	public static void display2DArray(double[][] Mat) {
		
		int rows = Mat.length;
		int cols = Mat[0].length;
		
		for (int i = 0; i < rows; i++) {
			
		    for (int j = 0; j < cols; j++) {
		    	
		        System.out.print(Mat[i][j] + " ");
		        
		    }
		    
		    System.out.println();
		}
	}
	
	public static void scalarByRealMatrix(RealMatrix aRealMat, double scalar) {
		int rows = aRealMat.getRowDimension();
		int cols = aRealMat.getColumnDimension();
		
		for (int i = 0; i < rows; i++) {
		    for (int j = 0; j < cols; j++) {
		    	aRealMat.multiplyEntry(i, j, scalar);
		    }
		}
	}
	
	public static void displayRealVector(RealVector vec) {
		
		int dim = vec.getDimension();
		
		for (int i = 0; i < dim; i++) { 
		        System.out.println(vec.getEntry(i)); 
		}
	}
	
	// For computing the EB covariance matrix and OU multivariate case
	public static void kronecker(double[][] A, double[][] B, double[][] Result) {
        int r1 = A.length;
        int c1 = A[0].length;
        int r2 = B.length;
        int c2 = B[0].length;

        //Result = new double[r1 * r2][c1 * c2]; Result matrix must have dimensions r1 * r2, c1 * c2
        for (int i1 = 0; i1 < r1; i1++) {
            for (int j1 = 0; j1 < c1; j1++) {
                for (int i2 = 0; i2 < r2; i2++) {
                    for (int j2 = 0; j2 < c2; j2++) {
                    	Result[i1 * r2 + i2][j1 * c2 + j2] = A[i1][j1] * B[i2][j2];
                    }
                }
            }
        }
    }
	
	public static void elementWiseProduct(double[][] A, double[][] B, double[][] Result) {
		
		int rows = A.length;
		int cols = A[0].length;
		
		for (int i = 0; i < rows; i++) {
			
		    for (int j = 0; j < cols; j++) {
		    	
		        Result[i][j] = A[i][j] * B[i][j];  
		    }
		}
	}
	
	public static int[] forLoopRange(int from, int limit) {
	    int[] numbers = new int[limit] ;
	    for (int i = 0; i < limit; i++) {
	        numbers[i] = from;
	        from++;
	    }
	    
	    return numbers;
	}

	
	public static void setBlocksInMatrix(double[][] blockMat, int m, int rowIndex, int colIndex, double[][] Result) {
		
		int i = 0, j = 0;
		int[] rowElements = forLoopRange(rowIndex * m, m); System.out.println(" Numero de filas filas " + rowElements.length);
		int[] colElements = forLoopRange(colIndex * m, m); System.out.println(" Numero de columnas " + colElements.length);
		
		for(int row:rowElements) {
			
			for(int col:colElements) {
				
				System.out.println(col);
				Result[row][col] = 	blockMat[i][j]; System.out.println("row = " + row + "," + "col = " + col);
				
				j++;
			}
	
			j = 0;
			i++;
		}
	}
	
	

}
