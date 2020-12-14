package contraband.utils;

import org.apache.commons.math3.linear.*;
import sun.security.krb5.Realm;

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
        SingularValueDecomposition xswsvd = new SingularValueDecomposition(xsw);
        double[] svdSingularValues = xswsvd.getSingularValues();
        RealMatrix leftSingular = xswsvd.getU();
        RealMatrix rightSingularTranspose = xswsvd.getVT();
        double max = Arrays.stream(svdSingularValues).max().getAsDouble();
        double tol = Math.max(n, p) * max * machineEps;
        int validNr = 0;
        for (int i = 0; i < svdSingularValues.length; i++) {
            if (svdSingularValues[i] > tol) {
                validNr += 1;
            }
        }
/*
        // start fast svd
        if (n > edgeRatio * p) {
            //psmall.svd
            xswsvd = pSmallSVD(xsw);
        } else if (edgeRatio * n < p) {
            //nsmall.svd
            xswsvd = nSmallSVD(xsw);
        } else{
            //positive.svd
            xswsvd = positiveSVD(xsw);
        }
*/

        // sweep(xswsvd$u, 2, xswsvd$d^3, "*")
        RealMatrix aMat = new Array2DRowRealMatrix(new double [n][validNr]);
        int k = 0;
        for (int j = 0; j < svdSingularValues.length; j ++) {
            if(svdSingularValues[j] > tol) {
                for (int i = 0; i < n; i++) {
                    aMat.setEntry(i, k, leftSingular.getEntry(i, j) * Math.pow(svdSingularValues[k], 3));
                }
                k ++;
            }
        }
        // (sweep(xswsvd$u, 2, xswsvd$d^3, "*") %*% t(xswsvd$v))
        RealMatrix dMat = new Array2DRowRealMatrix(new double[validNr][p]);
        k = 0;
        for (int j = 0; j < svdSingularValues.length; j ++) {
            if (svdSingularValues[j] > tol) {
                dMat.setRow(k, rightSingularTranspose.getRow(k));
                k ++;
            }
        }
        RealMatrix bMat = aMat.multiply(dMat);

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
                //sum1 += bMat.getEntry(i, j) * xsw.getEntry(i, j);
                sum3 += xsw.getEntry(i, j) * xsw.getEntry(i, j);
            }
            sum2 += sum3 * sum3;
        }
        double denominator = sum1 - sum2;

        // sweep(xs^2, MARGIN = 1, STATS = sw, FUN = "*")
        //RealMatrix xsSq = xs.preMultiply(xs);
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
                //sum4 += dMat.getEntry(i,j);
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

    public static SingularValueDecomposition positiveSVD(RealMatrix m){
        return new SingularValueDecomposition(m);
    }

    public static SingularValueDecomposition nSmallSVD (RealMatrix m){
        RealMatrix bMat = m.multiply(m.transpose());
        return new SingularValueDecomposition(bMat);

        //v = crossprod(m, u) %*% diag(1/d, nrow=length(d))
        // crossprod(m, u) = t(m) %*% u
    }

    public static SingularValueDecomposition pSmallSVD (RealMatrix m){
        RealMatrix bMat = m.transpose().multiply(m);
        return new SingularValueDecomposition(bMat);
    }

    public static List<Double> pSmallSingularValues (SingularValueDecomposition svd, int p){
        double[] svdSingularValues = svd.getSingularValues();
        double max = Arrays.stream(svdSingularValues).max().getAsDouble();
        double tol = p * max * machineEps;
        return getValidSingularValues(svdSingularValues, tol);
    }

    public static List<Double> getValidSingularValues(double[] singularValues, double tol){
        List<Double> valid = new ArrayList<>();
        int j = 0;
        for (double value : singularValues) {
            if (value > tol) {
                valid.add(j, value);
                j++;
            }
        }
        return valid;
    }




}
