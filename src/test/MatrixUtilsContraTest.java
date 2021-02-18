package test;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;

import contraband.math.MatrixUtilsContra;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Pau Bravo, Fabio K. Mendes, Rong Zhang
 */

public class MatrixUtilsContraTest {

	final static double EPSILON = 1e-4;

	@Test
	public void test1() {
		// block for matrix add
		// test a 1 * 1 matrix
		RealMatrix inRM1 = new Array2DRowRealMatrix(new double [][] {{ 2.0 }});
		RealMatrix rmToAdd1 = new Array2DRowRealMatrix(new double [][] {{3.0}});
		RealMatrix resRM1 = new Array2DRowRealMatrix(new double [1][1]);

		resRM1 = MatrixUtilsContra.matrixAdd(inRM1, rmToAdd1, resRM1);

		Assert.assertEquals(5, resRM1.getEntry(0,0), EPSILON);
	}

	@Test
	public void test2() {
		// test a 4 * 3 matrix
		RealMatrix inRM2 = new Array2DRowRealMatrix(new double [][]
				{{ 2.0, 3.0, 4.0 },
				 { 4.0, 5.0, 6.0 },
				 { -1.0, -2.0, -5.0 },
				 { 10.1, 9.2, -3.7 }
				});
		RealMatrix rmToAdd2 = new Array2DRowRealMatrix(new double [][]
				{{ -1.667465, -1.626329, -2.462770 },
				 { -2.329203, -1.973644, -1.833235 },
				 { -2.558182, -1.872955, -2.110709 },
				 { -2.274074, -2.245545, -2.314201 }
				});
		RealMatrix resRM2 = new Array2DRowRealMatrix(new double [4][3]);

		resRM2 = MatrixUtilsContra.matrixAdd(inRM2, rmToAdd2, resRM2);

		Assert.assertEquals(-6.014201, resRM2.getEntry(3,2), EPSILON);
	}

	@Test
	public void test3() {
		// block for matrix subtract
		RealMatrix inRM3 = new Array2DRowRealMatrix(new double [][]
				{{ -1.994546, -2.114393, -1.599438, -2.122436 },
				 { -1.735437, -1.670169, -1.780834, -2.245477 },
				 { -1.735442, -2.009275, -1.983139, -1.892316 },
				 { -1.692127, -1.942898, -1.601208, -1.982385 }
				});
		RealMatrix rmToSubtract1 = new Array2DRowRealMatrix(new double [][]
				{{ 3.930563, 2.338800, 2.682493, 3.316776 },
				 { 3.615184, 5.165212, 4.655430, 5.121321 },
				 { 3.125677, 2.930338, 3.339628, 4.388551 },
				 { 5.196331, 4.908884, 4.991581, 3.849233 }
				});
		RealMatrix resRM3 = new Array2DRowRealMatrix(new double [4][4]);

		resRM3 = MatrixUtilsContra.matrixSubtract(inRM3, rmToSubtract1, resRM3);

		Assert.assertEquals(-4.939613, resRM3.getEntry(2,1), EPSILON);
	}

	@Test
	public void test4() {
		// block for matrix scalar multiply
		RealMatrix inRM4 = new Array2DRowRealMatrix(new double [][]
				{{1.767012, 2.960950, 1.681112, 2.720202, 1.891744},
				 {1.319742, 2.048568, 1.982075, 3.104222, 1.966111},
				 {1.985207, 2.046110, 2.056098, 1.015987, 2.133164},
				 {1.879616, 2.094270, 2.605470, 2.279227, 2.642912}
				});
		double scalar1 = -2.0;
		RealMatrix resRM4 = new Array2DRowRealMatrix(new double [4][5]);

		resRM4 = MatrixUtilsContra.matrixScalarMultiply(inRM4, scalar1, resRM4);

		Assert.assertEquals(-5.285823, resRM4.getEntry(3,4), EPSILON);
	}

	@Test
	public void test5() {
		// block for matrix scalar add
		RealMatrix inRM5 = new Array2DRowRealMatrix(new double [][]
				{{2.550531, 2.350460, 1.242473, 1.534543},
				 {1.702243, 2.455985, 2.423319, 1.267112},
				 {2.077579, 2.054758, 2.861430, 2.958108},
				 {2.274022, 1.945004, 1.524684, 1.960117},
				 {1.887197, 1.782756, 1.224774, 1.812466},
		         {2.608558, 2.275520, 1.703467, 2.377170}
				});
		double scalar2 = 3.2;
		RealMatrix resRM5 = new Array2DRowRealMatrix(new double [6][4]);

		resRM5 = MatrixUtilsContra.matrixScalarAdd(inRM5, scalar2, resRM5);

		Assert.assertEquals(5.550460, resRM5.getEntry(0,1), EPSILON);
	}

	@Test
	public void test6() {
		// block for matrix transpose
		RealMatrix inRM6 = new Array2DRowRealMatrix(new double [][]
				{{5.555278, 3.637867, 3.826537, 5.580623, 4.233251, 4.265012, 5.668205, 5.960273},
				 {7.102671, 3.008732, 3.668539, 7.797644, 7.149419, 4.708736, 6.938855, 6.028421},
				 {8.244949, 2.804517, 5.590706, 4.507473, 3.285521, 5.566429, 1.413961, 6.178505},
				 {6.623658, 7.340270, 2.815117, 5.165536, 6.636394, 8.104003, 4.286587, 3.882813},
				 {5.732623, 1.996666, 6.853803, 5.431592, 1.740781, 7.074177, 6.466281, 2.540964},
				 {6.371464, 8.338275, 4.371867, 8.827669, 7.309662, 5.661179, 8.103298, 2.239197},
				 {1.422262, 1.983072, 4.297063, 4.046804, 5.934354, 7.656025, 5.650181, 2.041209}
				});
		RealMatrix resRM6 = new Array2DRowRealMatrix(new double [8][7]);

		resRM6 = MatrixUtilsContra.matrixTranspose(inRM6, resRM6);

		Assert.assertEquals(5.431592, resRM6.getEntry(3,4), EPSILON);
	}

	@Test
	public void test7() {
		RealMatrix inRM7 = new Array2DRowRealMatrix(new double [][]
				{{2.380636, 2.433278, 3.357671},
				 {3.346462, 2.153579, 3.418715},
				 {1.607831, 3.758575, 2.292255},
				 {2.131158, 3.110730, 2.237507}
				});
		RealMatrix rmToMultiply = new Array2DRowRealMatrix(new double [][]
				{{0.7944735, 0.3748316, 0.225552, 0.9200768, 0.9144262},
				 {2.6796039, 0.8591790, 2.216772, 3.0319734, 2.7087573},
				 {0.8331090, 2.0754329, 2.853639, 1.6177571, 1.1335097}
				});
		RealMatrix resRM7 = new Array2DRowRealMatrix(new double [4][5]);

		resRM7 = MatrixUtilsContra.matrixMultiply(inRM7, rmToMultiply, resRM7);

		Assert.assertEquals(15.23582, resRM7.getEntry(2,2), EPSILON);
	}


	/*
	 * In R:
	 *
	 * a <- matrix(c(1,2,3,4,5,6,7,8,9), nrow=3)
	 * b <- c(1,2,3)
	 *
	 * t(b)
	 *      [,1] [,2] [,3]
     * [1,]    1    2    3
	 *
	 * t(t(b))
	 *
	 *      [,1]
	 * [1,]    1
	 * [2,]    2
	 * [3,]    3
	 *
	 * t(b) %*% a %*% t(t(b))
	 *      [,1]
	 * [1,]  228
	 */
	@Test
	public void testTVecDotMatrixDotVec() {
		double[] inArr = new double[] { 1, 2, 3 };
		double[] inMatArr = new double[] { 1, 4, 7, 2, 5, 8, 3, 6, 9 };

		int nCol = 3;

		double res = MatrixUtilsContra.tVecDotMatrixDotVec(inArr, inMatArr, nCol);

		Assert.assertEquals(228, res, 0.0);
	}

	/*
	 * In R:
	 *
	 * a <- matrix(c(1,2,3,4,5,6,7,8,9), nrow=3)
	 * b <- c(1,2,3)
	 * t(b)
	 *       [,1] [,2] [,3]
	 * [1,]    1    2    3
	 *
	 * t(b) %*% a
	 * [,1] [,2] [,3]
	 * [1,]   14   32   50
	 */
	@Test
	public void testMatrixPreMultiply() {
		double[] inArr = new double[] { 1, 2, 3 };
		double[] inMatArr = new double[] { 1, 4, 7, 2, 5, 8, 3, 6, 9 };
		double[] resArr = new double[3];

		int nCol = 3;
		int nRow = 3;

		Double[] resArrDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixPreMultiply(inArr, inMatArr, nCol, nRow, resArr));

		Assert.assertArrayEquals(resArrDouble, new Double[] { 14.0, 32.0, 50.0 });
	}

	/*
	 * In R:
	 *
	 * b <- c(1,2,3)
	 * t(b)
	 *       [,1] [,2] [,3]
	 * [1,]    1    2    3
	 *
	 * 13 * t(b)
	 * [1] 13 26 39
	 */
	@Test
	public void testVectorPreMultiply() {
		double[] inArr = new double[] { 1, 2, 3 };
		double scalar = 13.0;
		double[] resArr = new double[3];

		Double[] resArrDouble = ArrayUtils.toObject(MatrixUtilsContra.vectorMapMultiply(inArr, scalar, resArr));

		Assert.assertArrayEquals(resArrDouble, new Double[] { 13.0, 26.0, 39.0 });
	}

	/*
	 * In R:
	 * a <- matrix(c(1.0, 5.0, 9.0, 13.0,
	 *               2.0, 6.0, 10.0, 14.0,
	 *               3.0, 7.0, 11.0, 15.0,
	 *               4.0, 8.0, 12.0, 16.0), ncol=4, byrow=T)
	 * b <- matrix(c(17.0, 21.0, 25.0, 29.0,
	 *               18.0, 22.0, 26.0, 30.0,
	 *               19.0, 23.0, 27.0, 31.0,
	 *               20.0, 24.0, 28.0, 32.0), ncol=4, byrow=T)
	 *
	 * a %*% b
	 *      [,1] [,2] [,3] [,4]
	 * [1,]  538  650  762  874
	 * [2,]  612  740  868  996
	 * [3,]  686  830  974 1118
	 * [4,]  760  920 1080 1240
	 */
	@Test
	public void testMatrixMultiply() {
		int nRowInMat = 4;
		int nColInMat = 4;
		int length = nRowInMat * nColInMat;

		double[] aMat = new double[] {
				1.0, 5.0, 9.0, 13.0,
				2.0, 6.0, 10.0, 14.0,
				3.0, 7.0, 11.0, 15.0,
				4.0, 8.0, 12.0, 16.0
		};

		double[] bMat = new double[] {
				17.0, 21.0, 25.0, 29.0,
				18.0, 22.0, 26.0, 30.0,
				19.0, 23.0, 27.0, 31.0,
				20.0, 24.0, 28.0, 32.0
		};

		double[] resMat = new double[length];

		Double[] resMatDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixMultiply(aMat, bMat, nRowInMat, nColInMat, resMat));

		Assert.assertArrayEquals(new Double[] {
				538.0, 650.0, 762.0, 874.0,
				612.0, 740.0, 868.0, 996.0,
				686.0, 830.0, 974.0, 1118.0,
				760.0, 920.0, 1080.0, 1240.0 }, resMatDouble);
	}

	/*
	 * In R
	 * a = matrix(c(3.92476468143533, 2.95493862847269, 5.47315682161247, 2.85610318223451, 3.62535769121946, 3.72081122250077, 3.70889115977579, 3.64833976126344, 4.32677785441579, 3.77038177701209, 3.31385684718126, 3.67462891649159, 3.73658224109536, 4.00675227771716, 3.23330588599873, 3.80698366724091), 4)
	 * a = t(a)
	 *
	 * res = t(a)
	 */
	@Test
	public void testVectorTranspose() {
		int ncol = 4;
		double[] aMat = new double[]{
				3.92476468143533, 2.95493862847269, 5.47315682161247, 2.85610318223451,
				3.62535769121946, 3.72081122250077, 3.70889115977579, 3.64833976126344,
				4.32677785441579, 3.77038177701209, 3.31385684718126, 3.67462891649159,
				3.73658224109536, 4.00675227771716, 3.23330588599873, 3.80698366724091
		};
		double[] resMat = new double[ncol * ncol];

		Double[] resMatDouble = ArrayUtils.toObject(MatrixUtilsContra.matrixTranspose(aMat, ncol, resMat));

		Assert.assertArrayEquals(new Double[]{3.92476468143533, 3.62535769121946, 4.32677785441579, 3.73658224109536, 2.95493862847269, 3.72081122250077, 3.77038177701209, 4.00675227771716, 5.47315682161247, 3.70889115977579, 3.31385684718126, 3.23330588599873, 2.85610318223451, 3.64833976126344, 3.67462891649159, 3.80698366724091}, resMatDouble);
	}

	/*
	 * In R
	 *
	 * a = c(2, -5, 3, -4)
     * b = c(3.95493862847269, 4.72081122250077, 4.77038177701209, 5.00675227771716)
	 * options(digits=17)
	 * res = a %*% b
	 */
	@Test
	public void testVectorDotMultiply(){
		double[] aVec = new double [] {2, -5, 3, -4};
		double[] bVec = new double [] {3.95493862847269, 4.72081122250077, 4.77038177701209, 5.00675227771716};

		double res = MatrixUtilsContra.vectorDotMultiply(aVec, bVec);

		Assert.assertEquals(-21.410042635390845, res, 0.0);
	}

	/*
	 * In R
	 *
	 * a = c(-0.206399781564555, 1.28555407499791, 1.0553689032024, 1.49801714307272, 1.20285942488946)
	 * options(digits=17)
	 * res = 0.2 * t(a) %*% a
	 */
	@Test
	public void testVecTransScalarMultiply () {
		double[] aVec = new double [] {
				-0.206399781564555, 1.28555407499791, 1.0553689032024, 1.49801714307272, 1.20285942488946
		};
		double scalar = 0.2;
		double res = MatrixUtilsContra.vecTransScalarMultiply(aVec, scalar, 5);
		Assert.assertEquals(1.2999959656810844, res, 0.0);
	}

	/*
	 * In R
	 *
	 * a = c(2.1, -5.2, 3.3, -4.4)
	 * b = c(3.95493862847269, 4.72081122250077, 4.77038177701209, 5.00675227771716)
	 *
	 * a + b
	 */
	@Test
	public void testVectorAdd () {
		int ncol = 4;

		double[] aVec = new double [] { 2.1, -5.2, 3.3, -4.4 };
		double[] bVec = new double [] { 3.95493862847269, 4.72081122250077, 4.77038177701209, 5.00675227771716 };

		double[] resVec = new double[ncol];

		Double[] resVecDouble = ArrayUtils.toObject(MatrixUtilsContra.vectorAdd(aVec, bVec, resVec));
		RealVector aRealVec = new ArrayRealVector(aVec);
		RealVector bRealVec = new ArrayRealVector(bVec);
		RealVector resRealVec = aRealVec.add(bRealVec);

		Assert.assertArrayEquals(new Double[] { resRealVec.getEntry(0), resRealVec.getEntry(1), resRealVec.getEntry(2), resRealVec.getEntry(3) }, resVecDouble);
	}

	@Test
	public void testGetMatrixEntry () {
		int ncol = 4;
		double[] aMat = new double[]{
				3.92476468143533, 2.95493862847269, 5.47315682161247, 2.85610318223451,
				3.62535769121946, 3.72081122250077, 3.70889115977579, 3.64833976126344,
				4.32677785441579, 3.77038177701209, 3.31385684718126, 3.67462891649159,
				3.73658224109536, 4.00675227771716, 3.23330588599873, 3.80698366724091
		};

		RealMatrix aMatRM = new Array2DRowRealMatrix(new double[][]{
				{3.92476468143533, 2.95493862847269, 5.47315682161247, 2.85610318223451},
				{3.62535769121946, 3.72081122250077, 3.70889115977579, 3.64833976126344},
				{4.32677785441579, 3.77038177701209, 3.31385684718126, 3.67462891649159},
				{3.73658224109536, 4.00675227771716, 3.23330588599873, 3.80698366724091}
		});

		Assert.assertEquals(aMatRM.getEntry(2 ,3), MatrixUtilsContra.getMatrixEntry(aMat, 2, 3, ncol), 0.0);
	}

	@Test
	public void testSetMatrixEntry () {
		int ncol = 4;
		double[] aMat = new double[]{
				3.92476468143533, 2.95493862847269, 5.47315682161247, 2.85610318223451,
				3.62535769121946, 3.72081122250077, 3.70889115977579, 3.64833976126344,
				4.32677785441579, 3.77038177701209, 3.31385684718126, 3.67462891649159,
				3.73658224109536, 4.00675227771716, 3.23330588599873, 3.80698366724091
		};

		RealMatrix aMatRM = new Array2DRowRealMatrix(new double[][]{
				{3.92476468143533, 2.95493862847269, 5.47315682161247, 2.85610318223451},
				{3.62535769121946, 3.72081122250077, 3.70889115977579, 3.64833976126344},
				{4.32677785441579, 3.77038177701209, 3.31385684718126, 3.67462891649159},
				{3.73658224109536, 4.00675227771716, 3.23330588599873, 3.80698366724091}
		});

		double value = 6.8;

		MatrixUtilsContra.setMatrixEntry(aMat, 2,3, value, ncol);
		aMatRM.setEntry(2, 3, value);

		Assert.assertEquals(aMatRM.getEntry(2 ,3), MatrixUtilsContra.getMatrixEntry(aMat, 2, 3, ncol), 0.0);
	}

	@Test
	public void testGetMatrixRow () {
		int ncol = 4;
		double[] aMat = new double[]{
				3.92476468143533, 2.95493862847269, 5.47315682161247, 2.85610318223451,
				3.62535769121946, 3.72081122250077, 3.70889115977579, 3.64833976126344,
				4.32677785441579, 3.77038177701209, 3.31385684718126, 3.67462891649159,
				3.73658224109536, 4.00675227771716, 3.23330588599873, 3.80698366724091
		};

		RealMatrix aMatRM = new Array2DRowRealMatrix(new double[][]{
				{3.92476468143533, 2.95493862847269, 5.47315682161247, 2.85610318223451},
				{3.62535769121946, 3.72081122250077, 3.70889115977579, 3.64833976126344},
				{4.32677785441579, 3.77038177701209, 3.31385684718126, 3.67462891649159},
				{3.73658224109536, 4.00675227771716, 3.23330588599873, 3.80698366724091}
		});

		double[] rowValues = new double [ncol];
		MatrixUtilsContra.getMatrixRow(aMat, 2, ncol, rowValues);

		Assert.assertArrayEquals(ArrayUtils.toObject(aMatRM.getRow(2)), ArrayUtils.toObject(rowValues));
	}

	@Test
	public void testSetMatrixRow () {
		int ncol = 4;
		double[] aMat = new double[]{
				3.92476468143533, 2.95493862847269, 5.47315682161247, 2.85610318223451,
				3.62535769121946, 3.72081122250077, 3.70889115977579, 3.64833976126344,
				4.32677785441579, 3.77038177701209, 3.31385684718126, 3.67462891649159,
				3.73658224109536, 4.00675227771716, 3.23330588599873, 3.80698366724091
		};

		RealMatrix aMatRM = new Array2DRowRealMatrix(new double[][]{
				{3.92476468143533, 2.95493862847269, 5.47315682161247, 2.85610318223451},
				{3.62535769121946, 3.72081122250077, 3.70889115977579, 3.64833976126344},
				{4.32677785441579, 3.77038177701209, 3.31385684718126, 3.67462891649159},
				{3.73658224109536, 4.00675227771716, 3.23330588599873, 3.80698366724091}
		});

		double[] rowValuesToSet = new double[] {5.2, 3.4, 1.0, -3.3};

		MatrixUtilsContra.setMatrixRow(aMat, rowValuesToSet, 2, ncol);

		aMatRM.setRow(2, rowValuesToSet);

		double[] rowValuesAfterSet = new double[ncol];

		MatrixUtilsContra.getMatrixRow(aMat, 2, ncol, rowValuesAfterSet);

		Assert.assertArrayEquals(ArrayUtils.toObject(aMatRM.getRow(2)), ArrayUtils.toObject(rowValuesAfterSet));
	}

	/*
	 * In R
	 *
	 * mat1 = matrix(c(-0.645193129729565, 0.880570797248417, 0.000937980841386414, -1.78306742105359, 1.8317655401661, 2.15197306144742, 3.68505700435419, -0.187213727070362, 3.33722112079121, -1.67453987788875, 4.63209934001223, -2.41553835887061, -1.8394774215567, -0.773758841460353, -2.20595620083541, -1.62315749365233), 4, 4, byrow = TRUE)
     * mat2 = matrix(c(3.99805567671961, 6.04887957663691, 1.16848962532753, 15.1709410155511, 3.59285822147303, 0.472026139032096, 11.476699867587, 2.5242582638748, 4.44900658426291, 3.804995006522, 1.06528397061267, 3.83606878754494, 9.12724650422943, 1.38426978606731, 1.22835160465911, 0.489186025224626), 4, 4, byrow = TRUE)
     * mat3 = matrix(c(4.28190062055364, 1.68441273737699, 2.42093200329691, -1.1208069534041, -2.09593287901953, -2.47336649568751, 3.97912367014214, 2.71986614447087, 3.35117470705882, -0.0328969140537083, -0.952275835443288, 2.56110355723649, 1.70456507708877, 4.74741416750476, 4.08583232900128, 2.53012317698449), 4, 4, byrow = TRUE)
     * resMat = mat1 %*% mat2 %*% mat3
	 */
	@Test
	public void testMatricesProduct () {
		int dim = 4;
		double[] mat1 = new double[]{-0.645193129729565, 0.880570797248417, 0.000937980841386414, -1.78306742105359, 1.8317655401661, 2.15197306144742, 3.68505700435419, -0.187213727070362, 3.33722112079121, -1.67453987788875, 4.63209934001223, -2.41553835887061, -1.8394774215567, -0.773758841460353, -2.20595620083541, -1.62315749365233};
		double[] mat2 = new double[]{3.99805567671961, 6.04887957663691, 1.16848962532753, 15.1709410155511, 3.59285822147303, 0.472026139032096, 11.476699867587, 2.5242582638748, 4.44900658426291, 3.804995006522, 1.06528397061267, 3.83606878754494, 9.12724650422943, 1.38426978606731, 1.22835160465911, 0.489186025224626};
		double[] mat3 = new double[]{4.28190062055364, 1.68441273737699, 2.42093200329691, -1.1208069534041, -2.09593287901953, -2.47336649568751, 3.97912367014214, 2.71986614447087, 3.35117470705882, -0.0328969140537083, -0.952275835443288, 2.56110355723649, 1.70456507708877, 4.74741416750476, 4.08583232900128, 2.53012317698449};
		double[] resMat = new double[dim * dim];
		MatrixUtilsContra.matricesProduct(mat1, mat2, mat3, dim, resMat);


        double[] exp = new double[]{-45.06403260135806, -51.97661534921737, -102.93874272355504, -1.6010562207429473, 256.0439385280864, 209.5276179684534, 338.941178953013, 234.7863513901713, 17.248822619312435, 226.0948195331119, 418.336093653783, 210.1761187428046, -220.6606489920099, -189.00845493027572, -317.4101985709398, -159.57540460518032};
		Assert.assertArrayEquals(ArrayUtils.toObject(resMat), ArrayUtils.toObject(exp));
	}

	/*
	 * In R
	 *
	 * mat1 = matrix(c(-1.32693286984645, -0.777384191935383, 0.970891128840399, -1.06308913385854, 0.696217525370674, 4.29829879893116, 3.09004870437948, -4.34566422008789, 1.47080558106243, -1.3191609278233, -0.924808043569947, 2.31426816310241, -0.41400931916902, 0.786774066268247, 0.613442564599179, 0.87506518006657), 4, 4, byrow = TRUE)
	 * mat2 = matrix(c(1.56312928421735, 0.61015987783989, 1.12751404910605, 1.93484731497927, 2.17207979625084, 1.96345202092547, 0.594608698949081, 0.908924838069615, 1.72934307609121, 1.87164572666989, 1.08674991988326, 1.46112581756923, 1.99173810018892, 2.36534894833742, 3.80466483365789, 2.32102529242019), 4, 4, byrow = TRUE)
	 * mat3 = matrix(c(3.63361590880609, 6.05236153562969, 12.5300851737367, 2.89721308508888, 0.855293145868927, 3.90544085276729, 1.11432315316051, 2.73023071000353, 2.83765158150345, 0.195223935879767, 5.02045130336142, 28.6674343769988, 11.2585298990117, 1.07210776535794, 0.142358511220664, 1.61179136950523), 4, 4, byrow = TRUE)
	 * resMat = mat1 %*% mat2 %*% t(mat1) + mat3
	 */
	@Test
	public void testMatricesTransProductAdd(){
		int dim = 4;
		double[] mat1 = new double[]{-1.32693286984645, -0.777384191935383, 0.970891128840399, -1.06308913385854, 0.696217525370674, 4.29829879893116, 3.09004870437948, -4.34566422008789, 1.47080558106243, -1.3191609278233, -0.924808043569947, 2.31426816310241, -0.41400931916902, 0.786774066268247, 0.613442564599179, 0.87506518006657};
		double[] mat2 = new double[]{1.56312928421735, 0.61015987783989, 1.12751404910605, 1.93484731497927, 2.17207979625084, 1.96345202092547, 0.594608698949081, 0.908924838069615, 1.72934307609121, 1.87164572666989, 1.08674991988326, 1.46112581756923, 1.99173810018892, 2.36534894833742, 3.80466483365789, 2.32102529242019};
		double[] mat3 = new double[]{3.63361590880609, 6.05236153562969, 12.5300851737367, 2.89721308508888, 0.855293145868927, 3.90544085276729, 1.11432315316051, 2.73023071000353, 2.83765158150345, 0.195223935879767, 5.02045130336142, 28.6674343769988, 11.2585298990117, 1.07210776535794, 0.142358511220664, 1.61179136950523};
		double[] resMat = new double[dim * dim];
		MatrixUtilsContra.matricesTransProductAdd(mat1, mat2, mat3, dim, resMat);

		double[] exp = new double[]{11.357975007004693, -6.414755369328954, 4.924296447425259, -4.568171863038602, -21.19022983803426, -1.3745230613746333, 14.17334201780606, -3.0883278807409185, 0.397800083631795, 12.88627486000583, 11.003160282178696, 39.54885315492658, 3.4830479675639943, 23.15240743364314, 2.7575692642906438, 8.49828141993841};
		Assert.assertArrayEquals(ArrayUtils.toObject(resMat), ArrayUtils.toObject(exp));
	}



}
