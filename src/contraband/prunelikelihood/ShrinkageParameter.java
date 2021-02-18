package contraband.prunelikelihood;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import contraband.math.MatrixUtilsContra;
import contraband.utils.ShrinkageUtils;
import org.apache.commons.math3.linear.*;
import java.util.ArrayList;
import java.util.List;

public class ShrinkageParameter extends CalculationNode {
    final public Input<MorphologicalData> traitInput = new Input<>("trait","Morphological data set.", Input.Validate.REQUIRED);
    final public Input<RealParameter> weightInput = new Input<>("weight", "Weight of each species.");

    private double delta;
    private double storedDelta;

    private static double machineEps = 2.220446E-16;
    private static int edgeRatio = 2;

    private int nTraits;
    private int nSpecies;

    private MorphologicalData morphData;
    private double[] traitValuesArr;

    private double[] singularValueArr;
    private double[] uMatArr;
    private double[] vMatArr;
    private double[] weight;

    private boolean equalWeight;
    private RealMatrix xsw;

    @Override
    public void initAndValidate() {
        morphData = traitInput.get();
        traitValuesArr = morphData .getMorphologicalData();
        nSpecies = morphData .getSpeciesNr();
        nTraits = morphData .getTotalTraitNr();

        weight = new double[nSpecies];
        if(weightInput.get() == null){
            // if weight is not specified, equal weight is used.
            for (int i = 0; i < nSpecies; i++) {
                weight[i] = 1.0 / nSpecies;
            }
            equalWeight = true;
        } else {
            weight = weightInput.get().getDoubleValues();
            equalWeight = false;
        }

        initIntermediateValues();

        estimateShrinkageParameter();
    }

    private void initIntermediateValues() {
        xsw = new Array2DRowRealMatrix(new double[nSpecies][nTraits]);
    }

    public double getDelta() { return delta; }

    public void estimateShrinkageParameter(){

        double[] xs = wtScale(traitValuesArr, weight, nSpecies, nTraits);

        double[] xsw = new double[nSpecies * nTraits];
        double[] sw = new double[nSpecies];
        double h1w2 = prepareMatrixForSVD(nSpecies, nTraits, weight, xs, sw, xsw);

        performSVD(nSpecies, nTraits, xsw);

        double denominator = calculateDenominator(nSpecies, nTraits, singularValueArr, uMatArr, vMatArr, xsw);

        double numerator = calculateNumerator(nSpecies, nTraits, xs, sw, denominator);

        if (denominator == 0) {
            delta = 1.0;
        } else {
            delta = Math.max(0, Math.min(1, numerator / denominator * h1w2));
        }

    }

    private void wtMoments(double[] data, double[] wVec, double [] m, double [] v, int n, int p) {
        double h1;
        if (equalWeight) {
            h1 = n / (n - 1.0);
        } else {
            double sum3 = 0;
            for (int i = 0; i < nSpecies; i++) {
                sum3 += wVec[i] * wVec[i];
            }
            h1 = 1 / (1 - sum3);
        }

        for (int j = 0; j < p; j ++) {
            double sum1 = 0;
            double sum2 = 0;
            for (int i = 0; i < n; i ++) {
                double d = MatrixUtilsContra.getMatrixEntry(data, i, j, nTraits);
                sum1 = sum1 + d * wVec[i];
                sum2 = sum2 + d * d * wVec[i];
            }
            m[j] = sum1;
            v[j] = h1 * (sum2 - sum1 * sum1);
            if (v[j] < machineEps) {
                v[j] = 0.0;
            }
        }
    }

    private double[] wtScale(double[] data, double[] wVec, int n, int p) {
        //compute column means and variances
        double[] m = new double[p]; // mean to return
        double[] v = new double[p]; // variance to return
        wtMoments(data, wVec, m, v, n, p);

        // each column in data subtracts each column in m
       // RealMatrix x = new Array2DRowRealMatrix(new double [n][p]);
        double[] x = new double[n * p];
        for (int j = 0; j < p; j ++) {
            for (int i = 0; i < n; i ++) {
                double value = MatrixUtilsContra.getMatrixEntry(data, i, j, p) - m[j];
                MatrixUtilsContra.setMatrixEntry(x, i, j, value, p);
            }
        }

        // each column in x is divided by square root of each column in m
        for (int j = 0; j < p; j ++) {
            double sc = Math.sqrt(v[j]);
            if (sc == 0.0) {
                //x.setColumn(j, new double [n]);
                for (int i = 0; i < n; i++) {
                    MatrixUtilsContra.setMatrixEntry(x, i, j, 0, p);
                }
            } else {
                for (int i = 0; i < n; i++) {
                    //x.setEntry(i, j, x.getEntry(i, j) / sc );
                    double value = MatrixUtilsContra.getMatrixEntry(x, i, j, p) / sc;
                    MatrixUtilsContra.setMatrixEntry(x, i, j, value, p);
                }
            }
        }
        return x;
    }

    private double prepareMatrixForSVD(int n, int p, double[] wVec, double[] xs, double[] sw, double[] xsw){
        double h1w2;
        if(equalWeight){
            h1w2 = 1.0 / (nSpecies - 1.0);
            for (int i = 0; i < n; i++) {
                sw[i] = Math.sqrt(wVec[i]);
            }
        } else {
            double w2 = 0.0;
            for (int i = 0; i < n; i++) {
                w2 += wVec[i] * wVec[i];
                sw[i] = Math.sqrt(wVec[i]);
            }
            h1w2 = w2 / (1 - w2);
        }

        // xsw = sweep(xs, MARGIN = 1, STATS = sw, FUN = "*")
        //RealMatrix xsw = new Array2DRowRealMatrix(new double[n][p]);
        for (int j = 0; j < p; j++) {
            for (int i = 0; i < n; i++) {
                //xsw.setEntry(i, j, xs.getEntry(i, j) * sw.getEntry(i));
                double value = MatrixUtilsContra.getMatrixEntry(xs, i, j, p) * sw[i];
                MatrixUtilsContra.setMatrixEntry(xsw, i, j, value, p);
            }
        }

        return h1w2;
    }

    private void performSVD(int n, int p, double[] xswArr){
        // For other cases, we need to implement fast.svd or get the delta from corpcor library.
        SingularValueDecomposition xswsvd;

        //RealMatrix xsw = new Array2DRowRealMatrix(new double[n][p]);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < p; j++){
                xsw.setEntry(i, j, MatrixUtilsContra.getMatrixEntry(xswArr, i, j, p));
            }
        }

        // start fast svd
        List<Double> singularValues = new ArrayList<>();
        List<Integer> validIndex = new ArrayList<>();
        if (n > edgeRatio * p) {
            //psmall.svd
            xswsvd = ShrinkageUtils.pSmallSVD(xsw);
            ShrinkageUtils.populateSingularValues(xswsvd, p, singularValues, validIndex);
            singularValueArr = ShrinkageUtils.populateSqrtSingularValuesArr(singularValues);
            double[] vTArr = ShrinkageUtils.getValidVMatrixArr(xswsvd.getV(), validIndex, p);
            uMatArr = ShrinkageUtils.getUMatArrForPSmall(xswArr, vTArr, singularValueArr, n);
            vMatArr = new double[p * singularValues.size()];
            MatrixUtilsContra.matrixTranspose(vTArr, p, singularValues.size(), vMatArr);
        } else if (edgeRatio * n < p) {
            //nsmall.svd
            xswsvd = ShrinkageUtils.nSmallSVD(xsw);
            ShrinkageUtils.populateSingularValues(xswsvd, n, singularValues, validIndex);
            singularValueArr = ShrinkageUtils.populateSqrtSingularValuesArr(singularValues);
            uMatArr = ShrinkageUtils.getValidUMatrixArr(xswsvd.getU(), validIndex, nSpecies);
            double[] vTArr = ShrinkageUtils.getVMatArrForNSmall(xswArr, uMatArr, singularValueArr, p, n);
            vMatArr = new double[p * singularValueArr.length];
            MatrixUtilsContra.matrixTranspose(vTArr, p, singularValueArr.length, vMatArr);
        } else{
            //positive.svd
            xswsvd = ShrinkageUtils.positiveSVD(xsw);
            ShrinkageUtils.populateSingularValues(xswsvd, Math.max(n, p), singularValues, validIndex);
            singularValueArr = ShrinkageUtils.populateSingularValuesArr(singularValues);
            uMatArr = ShrinkageUtils.getValidUMatrixArr(xswsvd.getU(), validIndex, n);
            vMatArr = ShrinkageUtils.getValidVMatrixArr(xswsvd.getVT(), validIndex, p);
        }
    }

    private double calculateDenominator(int n, int p, double[] singularValues, double[] uMatArr, double[] vMatArr, double[] xsw){
        int singularValueNr = singularValues.length;
        double[] aMat = new double[n * singularValueNr];
        for (int j = 0; j < singularValueNr; j ++) {
            for (int i = 0; i < n; i++) {
                double value = MatrixUtilsContra.getMatrixEntry(uMatArr, i, j, singularValueNr) * Math.pow(singularValues[j], 3);
                MatrixUtilsContra.setMatrixEntry(aMat, i, j, value, singularValueNr);
            }
        }

        double[] bMat = new double[n * p];
        MatrixUtilsContra.matrixMultiply(aMat, vMatArr, n, singularValueNr, p, bMat);

        double[] cMat = new double[n * p];
        for (int i = 0; i < n; i++) {
            for (int j =0; j < p; j++) {
                double value = MatrixUtilsContra.getMatrixEntry(xsw, i, j, p) * MatrixUtilsContra.getMatrixEntry(bMat, i, j, p);
                MatrixUtilsContra.setMatrixEntry(cMat, i, j, value, p);
            }
        }

        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int j = 0; j < p; j ++) {
            double sum3 = 0.0;
            for (int i = 0; i < n; i ++) {
                sum1 += MatrixUtilsContra.getMatrixEntry(cMat, i, j, p);
                sum3 += MatrixUtilsContra.getMatrixEntry(xsw, i, j, p) * MatrixUtilsContra.getMatrixEntry(xsw, i, j, p);
            }
            sum2 += sum3 * sum3;
        }
        return sum1 - sum2;
    }

    private double calculateNumerator(int n, int p, double[] xs, double[] sw, double denominator){
        double[] xs2w = new double[n * p];
        for (int j = 0; j < p; j ++) {
            for (int i = 0; i < n; i ++) {
                double xsSq = MatrixUtilsContra.getMatrixEntry(xs, i, j, p) * MatrixUtilsContra.getMatrixEntry(xs, i, j, p);
                MatrixUtilsContra.setMatrixEntry(xs2w, i, j, xsSq * sw[i], p);
            }
        }

        double[] xs2wPartial1 = new double[n * (p - 1)];
        for (int j = 0; j < (p -1); j ++) {
            for (int i = 0; i < n; i ++) {
                double value = MatrixUtilsContra.getMatrixEntry(xs2w, i, p - 2 - j, p);
                MatrixUtilsContra.setMatrixEntry(xs2wPartial1, i, j, value, p-1);
            }
        }

        double[] xs2wPartial2 = new double[n * (p - 1)];
        for (int j = 0; j < (p - 1); j ++) {
            for (int i = 0; i < n; i ++) {
                double value = MatrixUtilsContra.getMatrixEntry(xs2w, i, p - 1 - j, p);
                MatrixUtilsContra.setMatrixEntry(xs2wPartial2, i, j, value, p-1);
            }
        }

        for (int i = 0; i < n; i++) {
            double cumSum = 0.0;
            for (int j = 0; j < (p - 1); j ++) {
                cumSum += MatrixUtilsContra.getMatrixEntry(xs2wPartial2, i, j, p-1);
                MatrixUtilsContra.setMatrixEntry(xs2wPartial2, i, j, cumSum, p - 1);
            }
        }

        double sum4 = 0.0;
        for (int i = 0; i < n; i ++) {
            for (int j = 0; j < (p - 1); j++){
                sum4 += MatrixUtilsContra.getMatrixEntry(xs2wPartial1, i, j, p - 1) * MatrixUtilsContra.getMatrixEntry(xs2wPartial2, i, j, p - 1);
            }
        }
        double sER2 = 2 * sum4;

        return sER2 - denominator;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    public void store() {
        super.store();
        storedDelta = delta;
    }

    @Override
    public void restore() {
        super.restore();
        double tempDelta = delta;
        delta = storedDelta;
        storedDelta = tempDelta;

    }
}
