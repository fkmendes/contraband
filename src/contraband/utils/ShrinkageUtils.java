package contraband.utils;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/*
 * This class implements estimate.lambda() in R package corpcor
 * which returns shrinkage parameter given a matrix of continuous data set and a vector of weight
 */
public class ShrinkageUtils {

    private static double machineEps = 2.220446E-16;
    private static int edgeRatio = 2;
    /*
     * traitMat is an n * p RealMatrix
     * wVec is an n-dimensional vector
     */
    public static double estimateDelta(RealMatrix traitsMat, double[] wVec) {
        // number of species
        int n = traitsMat.getRowDimension();

        // number of traits
        int p = traitsMat.getColumnDimension();


        RealMatrix xs = wtScale(traitsMat, wVec, n, p);

        // w2 = sum(w * w)
        // sw = sqrt(w)
        RealVector sw = new ArrayRealVector(new double[n]);
        double w2 = 0.0;
        for (int i = 0; i < n; i++) {
            w2 = w2 + wVec[i] * wVec[i];
            sw.setEntry(i, Math.sqrt(wVec[i]));
        }
        double h1w2 = w2 / (1 - w2);

        // xsw = sweep(xs, MARGIN = 1, STATS = sw, FUN = "*")
        RealMatrix xsw = new Array2DRowRealMatrix(new double[n][p]);
        for (int j = 0; j < p; j++) {
            for (int i = 0; i < n; i++) {
                xsw.setEntry(i, j, xs.getEntry(i, j) * sw.getEntry(i));
            }
        }

        //xswsvd = fast.svd(xsw)
        // TO DO: this part currently only applies to the cases where nTraits and nSpecies are approximately the same.
        // For other cases, we need to implement fast.svd or get the delta from corpcor library.
        SingularValueDecomposition xswsvd;

        // start fast svd
        List<Double> singularValues = new ArrayList<>();
        List<Integer> validIndex = new ArrayList<>();
        RealMatrix vMat; RealMatrix uMat;
        if (n > edgeRatio * p) {
            //psmall.svd
            xswsvd = pSmallSVD(xsw);
            populateSingularValues(xswsvd, p, singularValues, validIndex);
            populateSqrtSingularValues(singularValues);
            vMat = getValidVMatrix(xswsvd.getVT(), validIndex);
            uMat = getUMatForPSmall(xsw, vMat, singularValues);
        } else if (edgeRatio * n < p) {
            //nsmall.svd
            xswsvd = nSmallSVD(xsw);
            populateSingularValues(xswsvd, n, singularValues, validIndex);
            populateSqrtSingularValues(singularValues);
            uMat = getValidUMatrix(xswsvd.getU(), validIndex);
            vMat = getVMatForNSmall(xsw, uMat, singularValues);
        } else{
            //positive.svd
            xswsvd = positiveSVD(xsw);
            populateSingularValues(xswsvd, Math.max(n, p), singularValues, validIndex);
            uMat = getValidUMatrix(xswsvd.getU(), validIndex);
            vMat = getValidVMatrix(xswsvd.getVT(), validIndex);
        }


        // sweep(xswsvd$u, 2, xswsvd$d^3, "*")
        RealMatrix aMat = new Array2DRowRealMatrix(new double [n][singularValues.size()]);
        for (int j = 0; j < singularValues.size(); j ++) {
                for (int i = 0; i < n; i++) {
                    aMat.setEntry(i, j, uMat.getEntry(i, j) * Math.pow(singularValues.get(j), 3));
                }
        }

        // (sweep(xswsvd$u, 2, xswsvd$d^3, "*") %*% t(xswsvd$v))
        RealMatrix bMat = aMat.multiply(vMat);

        // xsw * (sweep(xswsvd$u, 2, xswsvd$d^3, "*") %*% t(xswsvd$v))
        RealMatrix cMat = new Array2DRowRealMatrix(new double [n][p]);
        for (int i = 0; i < n; i++) {
            for (int j =0; j < p; j++) {
                cMat.setEntry(i, j, xsw.getEntry(i,j) * bMat.getEntry(i, j));
            }
        }

        //sum1 = sum(xsw * (sweep(xswsvd$u, 2, xswsvd$d^3, "*") %*% t(xswsvd$v)))
        double sum1 = 0.0;
        //sum2 = sum(colSums(xsw^2)^2)
        double sum2 = 0.0;
        for (int j = 0; j < p; j ++) {
            double sum3 = 0.0;
            for (int i = 0; i < n; i ++) {
                sum1 += cMat.getEntry(i, j);
                sum3 += xsw.getEntry(i, j) * xsw.getEntry(i, j);
            }
            sum2 += sum3 * sum3;
        }
        double denominator = sum1 - sum2;

        // sweep(xs^2, MARGIN = 1, STATS = sw, FUN = "*")
        RealMatrix xs2w = new Array2DRowRealMatrix(new double [n][p]);
        for (int j = 0; j < p; j ++) {
            for (int i = 0; i < n; i ++) {
                double xsSq =  xs.getEntry(i,j) * xs.getEntry(i,j);
                xs2w.setEntry(i,j, xsSq * sw.getEntry(i));
            }
        }

        //xs2wPartial1 = xs2w[, (p - 1):1]
        RealMatrix xs2wPartial1 = new Array2DRowRealMatrix(new double [n][p - 1]);
        for (int j = 0; j < (p -1); j ++) {
            xs2wPartial1.setColumn(j, xs2w.getColumn(p - 2 - j));
        }

        //xs2wPartial2 = xs2w[, p:2, drop = FALSE]
        //t(apply(xs2w[, p:2, drop = FALSE],1, cumsum))
        RealMatrix xs2wPartial2 = new Array2DRowRealMatrix(new double [n][p - 1]);
        for (int j = 0; j < (p - 1); j ++) {
            xs2wPartial2.setColumn(j, xs2w.getColumn(p - 1 - j));
        }
        for (int i = 0; i < n; i++) {
            double cumSum = 0.0;
            for (int j = 0; j < (p - 1); j ++) {
                cumSum += xs2wPartial2.getEntry(i, j);
                xs2wPartial2.setEntry(i, j, cumSum);
            }
        }
        
        //xs2wPartial2 = xs2wPartial2.transpose();
        double sum4 = 0.0;
        for (int i = 0; i < n; i ++) {
            for (int j = 0; j < (p - 1); j++){
                sum4 += xs2wPartial1.getEntry(i,j) * xs2wPartial2.getEntry(i,j);
            }
        }
        double sER2 = 2 * sum4;

        //double denominator = sE2R;
        //double numerator = sER2 - sE2R;
        double numerator = sER2 - denominator;

        if (denominator == 0) {
            return 1.0;
        } else {
            return Math.max(0, Math.min(1, numerator / denominator * h1w2));
        }

    }

    private static RealMatrix wtScale(RealMatrix data, double[] wVec, int n, int p) {

        //compute column means and variances
        double[] m = new double[p]; // mean to return
        double[] v = new double[p]; // variance to return
        wtMoments(data, wVec, m, v, n, p);

        // each column in data subtracts each column in m
        RealMatrix x = new Array2DRowRealMatrix(new double [n][p]);
        for (int j = 0; j < p; j ++) {
            for (int i = 0; i < n; i ++) {
                x.setEntry(i, j, data.getEntry(i, j) - m[j]);
            }
        }

        // each column in x is divided by square root of each column in m
        for (int j = 0; j < p; j ++) {
            double sc = Math.sqrt(v[j]);
            if (sc == 0.0) {
                x.setColumn(j, new double [n]);
            } else {
                for (int i = 0; i < n; i++) {
                    x.setEntry(i, j, x.getEntry(i, j) / sc );
                }
            }
        }

        return x;

    }

    private static void wtMoments(RealMatrix data, double[] wVec, double [] m, double [] v, int n, int p) {
        // if (uniform weights)
        double h1 = n /(n - 1.0);

        // else { h1 = 1/(1 - sum(wVec * wVec))}
        RealVector w = new ArrayRealVector(wVec);
        //RealMatrix wx = data.scalarMultiply(wVec[0]);
        //RealVector res = data.preMultiply(w);
       // RealMatrix wxSq = new Array2DRowRealMatrix(new double [n][p]);
        //for (int j = 0; j < p; j ++) {
            //for (int i = 0; i < n; i ++) {
                //wxSq.setEntry(i, j, data.getEntry(i, j) * data.getEntry(i, j) * wVec[j]);
            //}
       // }


        // sum of each column in wx
        //double[] m = new double[p]; // mean to return
        //double[] v = new double[p]; // variance to return
        for (int j = 0; j < p; j ++) {
            double sum1 = 0;
            double sum2 = 0;
            for (int i = 0; i < n; i ++) {
                sum1 = sum1 + data.getEntry(i, j) * wVec[i];
                sum2 = sum2 + data.getEntry(i, j) * data.getEntry(i, j) * wVec[i];
            }
            m[j] = sum1;
            v[j] = h1 * (sum2 - sum1 * sum1);
            if (v[j] < machineEps) {
                v[j] = 0.0;
            }
        }
    }

    // s = svd(m)
    public static SingularValueDecomposition positiveSVD(RealMatrix m){
        return new SingularValueDecomposition(m);
    }

    // B = m %*% t(m)     # nxn matrix
    // s = svd(B,nv=0)
    public static SingularValueDecomposition nSmallSVD (RealMatrix m){
        RealMatrix bMat = m.multiply(m.transpose());
        return new SingularValueDecomposition(bMat);
    }

    // B = crossprod(m)   # pxp matrix
    // s = svd(B,nu=0)
    public static SingularValueDecomposition pSmallSVD (RealMatrix m){
        RealMatrix bMat = m.transpose().multiply(m);
        return new SingularValueDecomposition(bMat);
    }


    //   tol = dim(B)[1]*max(s$d)*.Machine$double.eps
    //   Positive = s$d > tol
    //   For pSmall dim(B)[1] = p
    //   For nSmall dim(B)[1] = n
    public static void populateSingularValues(SingularValueDecomposition svd, int p, List<Double> values, List<Integer> index){
        double[] svdSingularValues = svd.getSingularValues();
        double max = Arrays.stream(svdSingularValues).max().getAsDouble();
        double tol = p * max * machineEps;
        populateValidSingularValues(svdSingularValues, tol, values, index);
    }

    public static void populateValidSingularValues(double[] singularValues, double tol, List<Double> values, List<Integer> index){
        int j = 0;
        for (int i = 0; i < singularValues.length; i++) {
            double value = singularValues[i];
            if (value > tol) {
                values.add(j, value);
                index.add(j, i);
                j++;
            }
        }
    }

    // For both pSmall and nSmall
    // d = sqrt(s$d[Positive])
    public static void populateSqrtSingularValues(List<Double> singularValues){
        for(int i = 0; i < singularValues.size(); i++){
            singularValues.set(i, FastMath.sqrt(singularValues.get(i)));
        }
    }

    // For pSmall v = s$v[, Positive, drop=FALSE]
    // For nSmall u = s$u[, Positive, drop=FALSE]
    public static RealMatrix getValidUMatrix(RealMatrix aMat, List<Integer> index){
        int nRow = aMat.getRowDimension();
        int nCol = index.size();
        RealMatrix resMat = new Array2DRowRealMatrix(new double[nRow][nCol]);
        for(int i = 0; i < nCol; i ++){
            resMat.setColumn(i, aMat.getColumn(index.get(i)));
        }
        return resMat;
    }

    public static RealMatrix getValidVMatrix(RealMatrix aMat, List<Integer> index){
        int nCol = aMat.getColumnDimension();
        int nRow = index.size();
        RealMatrix resMat = new Array2DRowRealMatrix(new double[nRow][nCol]);
        for(int i = 0; i < nRow; i ++){
            resMat.setRow(i, aMat.getRow(index.get(i)));
        }
        return resMat;
    }

    // u = m %*% v %*% diag(1/d, nrow=length(d))
    public static RealMatrix getUMatForPSmall(RealMatrix m, RealMatrix v, List<Double> d){
        RealMatrix diag = MatrixUtils.createRealIdentityMatrix(d.size());
        for(int i = 0; i < d.size(); i++){
            diag.setEntry(i, i, 1/d.get(i));
        }
        return m.multiply(v).multiply(diag);
    }

    //  v = crossprod(m, u) %*% diag(1/d, nrow=length(d))
    public static RealMatrix getVMatForNSmall(RealMatrix m, RealMatrix u, List<Double> d){
        RealMatrix diag = MatrixUtils.createRealIdentityMatrix(d.size());
        for(int i = 0; i < d.size(); i++){
            diag.setEntry(i, i, 1/d.get(i));
        }
        return m.transpose().multiply(u).multiply(diag);
    }




}
