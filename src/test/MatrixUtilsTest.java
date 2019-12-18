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
	private RealMatrix resRM2;
	private RealMatrix resRM3;
	private RealMatrix resRM4;
	private RealMatrix resRM5;
	private RealMatrix resRM6;
	@Before
	public void setUP() throws Exception {
		// block for matrix add
		// test a 1 * 1 matrix
		RealMatrix inRM1 = new Array2DRowRealMatrix(new double [][]
				{{2.0}}
		);
		RealMatrix rmToAdd1 = new Array2DRowRealMatrix(new double [][]
				{{3.0}}
		);
        System.out.println("Display [2.0] + [3.0] = ");
		resRM1 = new Array2DRowRealMatrix(new double [1][1]);
        resRM1 = MatrixUtils.matrixAdd(inRM1, rmToAdd1, resRM1);
        GeneralUtils.displayRealMatrix(resRM1);

        // test a 4 * 3 matrix
		System.out.println("Display matrix add resRM2:");
		RealMatrix inRM2 = new Array2DRowRealMatrix(new double [][]
				{{2.0, 3.0, 4.0},
				 {4.0, 5.0, 6.0},
				 {-1.0, -2.0, -5.0},
				 {10.1, 9.2, -3.7}
				});
		RealMatrix rmToAdd2 = new Array2DRowRealMatrix(new double [][]
				{{-1.667465, -1.626329, -2.462770},
				 {-2.329203, -1.973644, -1.833235},
				 {-2.558182, -1.872955, -2.110709},
				 {-2.274074, -2.245545, -2.314201}
				});
		resRM2 = new Array2DRowRealMatrix(new double [4][3]);
		resRM2 = MatrixUtils.matrixAdd(inRM2, rmToAdd2, resRM2);
		GeneralUtils.displayRealMatrix(resRM2);

		// block for matrix subtract
		System.out.println("Display matrix subtract resRM3:");
		RealMatrix inRM3 = new Array2DRowRealMatrix(new double [][]
				{{-1.994546, -2.114393, -1.599438, -2.122436,},
				 {-1.735437, -1.670169, -1.780834, -2.245477},
				 {-1.735442, -2.009275, -1.983139, -1.892316},
				 { -1.692127, -1.942898, -1.601208, -1.982385}
				});
		RealMatrix rmToSubtract1 = new Array2DRowRealMatrix(new double [][]
				{{3.930563, 2.338800, 2.682493, 3.316776},
				 {3.615184, 5.165212, 4.655430, 5.121321},
				 {3.125677, 2.930338, 3.339628, 4.388551},
				 { 5.196331, 4.908884, 4.991581, 3.849233}
				});
		resRM3 = new Array2DRowRealMatrix(new double [4][4]);
		resRM3 = MatrixUtils.matrixSubtract(inRM3, rmToSubtract1, resRM3);
		GeneralUtils.displayRealMatrix(resRM3);

		// block for matrix scalar multiply
		System.out.println("Display matrix scalar multiply resRM4:");
		RealMatrix inRM4 = new Array2DRowRealMatrix(new double [][]
				{{1.767012, 2.960950, 1.681112, 2.720202, 1.891744},
				 {1.319742, 2.048568, 1.982075, 3.104222, 1.966111},
				 {1.985207, 2.046110, 2.056098, 1.015987, 2.133164},
				 {1.879616, 2.094270, 2.605470, 2.279227, 2.642912}
				});
		double scalar1 = -2.0;
		resRM4 = new Array2DRowRealMatrix(new double [4][5]);
		resRM4 = MatrixUtils.matrixScalarMultiply(inRM4, scalar1, resRM4);
		GeneralUtils.displayRealMatrix(resRM4);

		// block for matrix scalar add
		System.out.println("Display matrix scalar add resRM5:");
		RealMatrix inRM5 = new Array2DRowRealMatrix(new double [][]
				{{2.550531, 2.350460, 1.242473, 1.534543},
				 {1.702243, 2.455985, 2.423319, 1.267112},
				 {2.077579, 2.054758, 2.861430, 2.958108},
				 {2.274022, 1.945004, 1.524684, 1.960117},
				 {1.887197, 1.782756, 1.224774, 1.812466},
		         {2.608558, 2.275520, 1.703467, 2.377170}
				});
		double scalar2 = 3.2;
		resRM5 = new Array2DRowRealMatrix(new double [6][4]);
		resRM5 = MatrixUtils.matrixScalarAdd(inRM5, scalar2, resRM5);
		GeneralUtils.displayRealMatrix(resRM5);

		// block for matrix transpose
		System.out.println("Display matrix transpose resRM6:");
		RealMatrix inRM6 = new Array2DRowRealMatrix(new double [][]
				{{5.555278, 3.637867, 3.826537, 5.580623, 4.233251, 4.265012, 5.668205, 5.960273},
				 {7.102671, 3.008732, 3.668539, 7.797644, 7.149419, 4.708736, 6.938855, 6.028421},
				 {8.244949, 2.804517, 5.590706, 4.507473, 3.285521, 5.566429, 1.413961, 6.178505},
				 {6.623658, 7.340270, 2.815117, 5.165536, 6.636394, 8.104003, 4.286587, 3.882813},
				 {5.732623, 1.996666, 6.853803, 5.431592, 1.740781, 7.074177, 6.466281, 2.540964},
				 {6.371464, 8.338275, 4.371867, 8.827669, 7.309662, 5.661179, 8.103298, 2.239197},
				 {1.422262, 1.983072, 4.297063, 4.046804, 5.934354, 7.656025, 5.650181, 2.041209}
				});
		resRM6 = new Array2DRowRealMatrix(new double [8][7]);
		resRM6 = MatrixUtils.matrixTranspose(inRM6, resRM6);
		GeneralUtils.displayRealMatrix(resRM6);

	}

	@Test
	public void againstresLikelihood () {
		Assert.assertEquals(5, resRM1.getEntry(0,0), EPSILON);
		Assert.assertEquals(-6.014201, resRM2.getEntry(3,2), EPSILON);
		Assert.assertEquals(-4.939613, resRM3.getEntry(2,1), EPSILON);
		Assert.assertEquals(-5.285823, resRM4.getEntry(3,4), EPSILON);
		Assert.assertEquals(5.550460, resRM5.getEntry(0,1), EPSILON);
		Assert.assertEquals(5.431592, resRM6.getEntry(3,4), EPSILON);
	}
}
