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
import contraband.clock.ColorManager;
import contraband.math.MVNUtils;

public class MVNShiftVcvMatTestDriver {

	public static void main(String[] args) {
		// tree
		String treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// VCV Mat
		RealParameter colorValues = new RealParameter(new Double[] { 0.1160941, 0.01914707 });
		IntegerParameter colorAssignments = new IntegerParameter(new Integer[] { 1, 1, 0, 0, 0 });
		
		ColorManager colors = new ColorManager();
		colors.initByName("nTraits", 1, "maxNColors", 2, "tree", myTree, "colorValues", colorValues, "colorAssignments", colorAssignments, "coalCorrection", false);
		
		double[][] colorValuesMat = colors.getSpColorValuesMatOneTrait();
		for (int i=0; i<colorValuesMat.length; ++i) {
			System.out.println(Arrays.toString(colorValuesMat[i]));
		}
		
		RealMatrix vcvMat = new Array2DRowRealMatrix(colorValuesMat);
		LUDecomposition vcvMatLUD = new LUDecomposition(vcvMat);
		double detVcvMat = vcvMatLUD.getDeterminant();
		RealMatrix invVcvMat = vcvMatLUD.getSolver().getInverse();
		
		// mean vector (root value repeated)
		double[] meanInput = new double[] { -1.125558, -1.125558, -1.125558 };
		RealVector mean = new ArrayRealVector(meanInput);
		
		// data
		RealVector dataVec = new ArrayRealVector(new double[] { -0.9291812, -0.7312343, -1.6712572 });
		
		double resLogLk = MVNUtils.getMVNLogLk(3, mean, dataVec, invVcvMat, detVcvMat);
		System.out.println(resLogLk); // -0.8583676
	}

}
