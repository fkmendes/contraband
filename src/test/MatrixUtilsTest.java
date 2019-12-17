package test;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import contraband.MatrixUtils;
import contraband.GeneralUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;


public class MatrixUtilsTest {
	final static double EPSILON = 1e-4;
	private RealMatrix resRM1;

	@Before
	public void setUP() throws Exception {
		RealMatrix inRM1 = new Array2DRowRealMatrix(new double [][]
				{{2.0}}
		);
		RealMatrix rmToAdd = new Array2DRowRealMatrix(new double [][]
				{{3.0}}
		);

		resRM1 = new Array2DRowRealMatrix(new double [1][1]);
        resRM1 = MatrixUtils.matrixAdd(inRM1, rmToAdd, resRM1);
        System.out.println(resRM1);
	}

	@Test
	public void againstresLikelihood () {
		Assert.assertEquals(5, resRM1.getEntry(0,0), EPSILON);
	}
}
