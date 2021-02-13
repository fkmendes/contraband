package testdrivers;

import com.sun.deploy.util.GeneralUtil;
import contraband.math.MatrixUtilsContra;
import contraband.prunelikelihood.OUPruneUtils;
import contraband.utils.GeneralUtils;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;

public class OUModelParamsTestDriver1 {
    public static void main(String[] args) {

        // (1) initialize the test
        int nTraits = 4;
        double[] res = new double[nTraits * nTraits];
        double[] phiTranspose = new double[nTraits * nTraits];
        double[] eMulPhiMat = new double[nTraits * nTraits];
        double[] resVec = new double[nTraits];
        double[] eTranspose = new double[nTraits * nTraits];

        // (2) specify the input: parameters from NodeMath
        double varianceMatrixDet = 10.828800180056;
        double[] invVariance = new double[] { 0.303159451068479, -0.177162352807323, -0.167490282668518, 0.463804357145631, -0.177162352807323, 0.145356310238897, 0.0876205939159736, -0.270101488601486, -0.167490282668518, 0.0876205939159733, 0.139138516499997, -0.17669132609597, 0.463804357145632, -0.270101488601486, -0.176691326095971, 0.889433801166107 };
        double[] phiMat = new double[] {5.79139780244504, -0.655933465833111, 4.23336767300133, -2.23847639400205, -0.655933465833112, 2.52293429844627, -2.23165633273757, 0.372164773841502, 4.23336767300133, -2.23165633273757, 5.6575229085707, -1.6836916603114, -2.23847639400205, 0.372164773841502, -1.6836916603114, 1.08811139295935};
        double[] omegaVec = new double[]{-18.6505144809713, 4.67709793742575, -16.945172311599, 6.30278768018728};
        double [] traitsVec = new double[] {-1107.45586477186, 348.317941028756, -1096.18225437494, 442.671195094954};

        RealMatrix invVarRM = new Array2DRowRealMatrix(new double[][]{
                {0.303159451068479, -0.177162352807323, -0.167490282668518, 0.463804357145631},
                {-0.177162352807323, 0.145356310238897, 0.0876205939159736, -0.270101488601486},
                {-0.167490282668518, 0.0876205939159733, 0.139138516499997, -0.17669132609597},
                {0.463804357145632, -0.270101488601486, -0.176691326095971, 0.889433801166107}
        });
        RealMatrix phiRM = new Array2DRowRealMatrix(new double[][]{
                {5.79139780244504, -0.655933465833111, 4.23336767300133, -2.23847639400205},
                {-0.655933465833112, 2.52293429844627, -2.23165633273757, 0.372164773841502},
                {4.23336767300133, -2.23165633273757, 5.6575229085707, -1.6836916603114},
                {-2.23847639400205, 0.372164773841502, -1.6836916603114, 1.08811139295935}
        });
        RealVector omegaRV = new ArrayRealVector(omegaVec);
        RealVector traitsRV = new ArrayRealVector(traitsVec);

        // (3) AbCdEf
        // aMat
        double[] aMat = new double [nTraits * nTraits];
        MatrixUtilsContra.vectorMapMultiply(invVariance, -0.5, aMat);
        RealMatrix aRM = OUPruneUtils.getAMatForOU(invVarRM);
        System.out.println("aMat = " + Arrays.toString(aMat) + "\n" + "A Matrix = ");
        GeneralUtils.displayRealMatrix(aRM);

        // eMat
        double[] eMat = new double [nTraits * nTraits];
        MatrixUtilsContra.matrixTranspose(phiMat, nTraits, phiTranspose);
        MatrixUtilsContra.matrixMultiply(phiTranspose, invVariance, nTraits, nTraits, eMat);
        RealMatrix eRM = OUPruneUtils.getEMatForOU(phiRM, invVarRM);
        System.out.println("eMat = " + Arrays.toString(eMat) + "\n" + "E Matrix = ");
        GeneralUtils.displayRealMatrix(eRM);

        // CMat
        double[] cMat = new double[nTraits * nTraits];
        MatrixUtilsContra.matrixMultiply(eMat, phiMat, nTraits, nTraits, eMulPhiMat);
        MatrixUtilsContra.vectorMapMultiply(eMulPhiMat, -0.5, cMat);
        RealMatrix cRM = OUPruneUtils.getCMatForOU(phiRM, eRM);
        System.out.println("cMat = " + Arrays.toString(cMat) + "\n" + "C Matrix = ");
        GeneralUtils.displayRealMatrix(cRM);

        // f
        double f = -0.5 * (MatrixUtilsContra.tVecDotMatrixDotVec(omegaVec, invVariance, nTraits) + nTraits * FastMath.log(2 * Math.PI) +  varianceMatrixDet);
        double f1 = OUPruneUtils.getFforOU(omegaRV, invVarRM, FastMath.exp(varianceMatrixDet), nTraits);
        System.out.println("f = " + f + ", f1 = "+ f1);

        // bVec
        double[] bVec = new double[nTraits];
        MatrixUtilsContra.matrixPreMultiply(omegaVec, invVariance, nTraits, nTraits, bVec);
        RealVector bRV = OUPruneUtils.getBVecForOU(invVarRM, omegaRV);
        System.out.println("bVec = " + Arrays.toString(bVec) + "\n" + "b vector = ");
        GeneralUtils.displayRealVector(bRV);

        // dVec
        double[] dVec = new double[nTraits];
        MatrixUtilsContra.matrixTranspose(eMat, nTraits, eTranspose);
        MatrixUtilsContra.matrixPreMultiply(omegaVec, eTranspose, nTraits, nTraits, resVec);
        MatrixUtilsContra.vectorMapMultiply(resVec, -1, dVec);
        RealVector dRV = OUPruneUtils.getDVecForOU(eRM, omegaRV);
        System.out.println("dVec = " + Arrays.toString(dVec) + "\n" + "d vector = ");
        GeneralUtils.displayRealVector(dRV);

        // (4) Lmr
        //L
        double[] lMat = new double [nTraits * nTraits];
        MatrixUtilsContra.matrixTransAddScalar(cMat, 0.5, nTraits, lMat);
        System.out.println("lMat = " + Arrays.toString(lMat) + "\n" + "L Matrix = ");
        RealMatrix lRM = OUPruneUtils.getLMatForOULeaf(cRM);
        GeneralUtils.displayRealMatrix(lRM);

        //r
        double r = MatrixUtilsContra.tVecDotMatrixDotVec(traitsVec, aMat, nTraits)
                + MatrixUtilsContra.vectorDotMultiply(traitsVec, bVec)
                + f;
        double r1 = OUPruneUtils.getRForOULeaf(aRM, traitsRV, bRV, f);
        System.out.println("r = " + r + ", r1 = "+ r1);

        // m
        double[] mVec = new double[nTraits];
        MatrixUtilsContra.matrixTranspose(eMat, nTraits, eTranspose);
        MatrixUtilsContra.matrixPreMultiplyAddVector(traitsVec, eTranspose, dVec, nTraits, nTraits, mVec);
        RealVector mRV = OUPruneUtils.getMVecForOULeaf(eRM, traitsRV, dRV);
        System.out.println("mVec = " + Arrays.toString(mVec) + "\n" + "m vector = ");
        GeneralUtils.displayRealVector(mRV);
    }
}
