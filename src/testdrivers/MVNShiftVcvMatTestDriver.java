package testdrivers;

import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.ColorManager;
import contraband.MVNUtils;

public class MVNShiftVcvMatTestDriver {

	public static void main(String[] args) {
		// tree
		String treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.2804859, 0.226648 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 1, 1, 0, 0, 0 });
		
		ColorManager colors = new ColorManager();
		colors.initByName("tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments);
		
		double[][] colorValuesMat = colors.getSpColorValuesMatrix();
		for (int i=0; i<colorValuesMat.length; ++i) {
			System.out.println(Arrays.toString(colorValuesMat[i]));
		}
		
		RealMatrix vcvMat = new Array2DRowRealMatrix(colorValuesMat);
		LUDecomposition vcvMatLUD = new LUDecomposition(vcvMat);
		double detVcvMat = vcvMatLUD.getDeterminant();
		RealMatrix invVcvMat = vcvMatLUD.getSolver().getInverse();
		
		// mean vector (root value repeated)
		double[] meanInput = new double[] { 0.362281, 0.362281, 0.362281 };
		RealVector mean = new ArrayRealVector(meanInput);
		
		// data
		RealVector dataVec = new ArrayRealVector(new double[] { 0.1466427, -0.5456044, 1.1624960 });
		
		double resLogLk = MVNUtils.getMVNLogLk(3, mean, dataVec, invVcvMat, detVcvMat);
		System.out.println(resLogLk); // -3.106220692988152
	}

}
