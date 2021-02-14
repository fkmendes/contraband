package testdrivers;

import contraband.math.MatrixUtilsContra;
import contraband.prunelikelihood.OUPruneUtils;
import contraband.utils.GeneralUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import java.util.Arrays;

public class OUModelParamsTestDriver2 {
    public static void main(String[] args) {

        // (1) initialize the test
        int nTraits = 4;
        double[] bPlusM = new double[nTraits];
        double[] dToSubtract = new double[nTraits];

        // (2) specify the input: parameters from NodeMath
        double[] aPlusLInverse = new double[] {-7.60493184544515, -6.03637968255101, -1.15913549831401, 1.09530716933721, -6.03637968255101, -12.9469072534224, -6.14516461715909, -2.52043148954421, -1.15913549831401, -6.14516461715909, -5.54737311159581, -2.71427244168607, 1.09530716933721, -2.5204314895442, -2.71427244168607, -3.88102821721889};
        RealMatrix aPlusLInvRM = new Array2DRowRealMatrix(new double[][]{
                {-7.60493184544515, -6.03637968255101, -1.15913549831401, 1.09530716933721},
                {-6.03637968255101, -12.9469072534224, -6.14516461715909, -2.5204314895442},
                {-1.15913549831401, -6.14516461715909, -5.54737311159581, -2.71427244168607},
                {1.09530716933721, -2.52043148954421, -2.71427244168607, -3.88102821721889}
        });
        RealMatrix eAPlusLInvRM = new Array2DRowRealMatrix(new double[][] {
                {0.118039645101024, 1.08598015604285, 0.652426637069064, 0.138301664296393},
                {-0.60446173227927, -0.88437136086398, -0.0749117251358719, 0.553576490289782},
                {0.992993122313137, 0.623549171651613, -0.367448119133317, -0.499030974060303},
                {-0.215510097837987, -0.57621692848297, -0.343781890606824, -0.416785481660141}
        });
        double logDetNegativeTwoAplusL = -2.7849175386693;
        double[] bVec = new double[]{-0.72127331852945, 0.796876328617564, 0.0622159908252355, -1.31350362352276};
        RealVector bRV = new ArrayRealVector(bVec);
        double[] mVec = new double[]{-85.2785852896564, 107.806637958107, -152.646008693384, 33.0371416900153};
        RealVector mRV = new ArrayRealVector(mVec);
        double rC = -13324.4944559487;
        double fC = -13.0132502046673;
        double[] eMat = new double[]{0.124660641703936, -0.145815733087733, -0.0429322745118151, 0.124268224350433, -0.0994290120574335, 0.186869823978932, -0.0453441556462196, -0.260342862825273, -0.0497327328190371, -0.12389556635717, 0.180086055278854, 0.0690605711006969, 0.0421239354253536, 0.00924367274113691, -0.0189938623244606, 0.126559400203708};
        RealMatrix eRM = new Array2DRowRealMatrix(new double[][]{
                {0.124660641703936, -0.145815733087733, -0.0429322745118151, 0.124268224350433},
                {-0.0994290120574335, 0.186869823978932, -0.0453441556462196, -0.260342862825273},
                {-0.0497327328190371, -0.12389556635717, 0.180086055278854, 0.0690605711006969},
                {0.0421239354253536, 0.00924367274113691, -0.0189938623244606, 0.126559400203708}
        });
        double[] dVec = new double[]{1.49624854496904, -1.85588943981415, 2.26824576514175, -0.377131793432348};
        RealVector dRV = new ArrayRealVector(dVec);
        double[] cMat = new double[]{-0.17884259790653, 0.153797131560893, -0.20051261764744, 0.0629074998951295, 0.153797131560893, -0.270490617090084, 0.328073641930157, -0.0425897018690009, -0.20051261764744, 0.328073641930157, -0.484259828744313, 0.0814239061005466, 0.0629074998951296, -0.0425897018690009, 0.0814239061005468, -0.039418633571801};
        RealMatrix cRM = new Array2DRowRealMatrix(new double[][]{
                {-0.17884259790653, 0.153797131560893, -0.20051261764744, 0.0629074998951295},
                {0.153797131560893, -0.270490617090084, 0.328073641930157, -0.042589701869000},
                {-0.20051261764744, 0.328073641930157, -0.484259828744313, 0.0814239061005466},
                {0.0629074998951296, -0.0425897018690009, 0.0814239061005468, -0.039418633571801}
        });

        //r
        MatrixUtilsContra.vectorAdd(bVec, mVec, bPlusM);
        double r = rC+ fC
                + 0.5 * nTraits* FastMath.log(2 * Math.PI)
                - 0.5 * logDetNegativeTwoAplusL
                - 0.25 * MatrixUtilsContra.tVecDotMatrixDotVec(bPlusM, aPlusLInverse, nTraits);
        double r1 = OUPruneUtils.getRForOUIntNode(bRV, mRV, aPlusLInvRM, fC, rC, nTraits, logDetNegativeTwoAplusL);
        System.out.println("r = " + r + ", r1 = " + r1);

        //m
        double[] eAPlusLInv = new double[nTraits * nTraits];
        MatrixUtilsContra.matrixMultiply(eMat, aPlusLInverse, nTraits, nTraits, eAPlusLInv);
        MatrixUtilsContra.matrixTransPreMapMultiply(bPlusM, eAPlusLInv, 0.5, nTraits, nTraits, dToSubtract);
        MatrixUtilsContra.vectorSubtract(dVec, dToSubtract, mVec);
        RealVector m = OUPruneUtils.getMVecForOUIntNode(eAPlusLInvRM, bRV, mRV, dRV);
        System.out.println("mVec = " + Arrays.toString(mVec) + "\n" + "m vector = ");
        GeneralUtils.displayRealVector(m);

        //L
        double[] eMatTransScalar = new double [nTraits * nTraits];
        double[] cToSubtract = new double [nTraits * nTraits];
        double[] lMat = new double [nTraits * nTraits];
        MatrixUtilsContra.matrixTransScalar(eMat, 0.25, nTraits, eMatTransScalar);
        MatrixUtilsContra.matrixMultiply(eAPlusLInv, eMatTransScalar, nTraits, nTraits, cToSubtract);
        MatrixUtilsContra.matrixSubtract(cMat, cToSubtract, nTraits, lMat);
        RealMatrix lRM = OUPruneUtils.getLMatForOUIntNode(cRM, eRM, eAPlusLInvRM);
        System.out.println("lMat = " + Arrays.toString(lMat) + "\n" + "L Matrix = ");
        GeneralUtils.displayRealMatrix(lRM);
    }
}
