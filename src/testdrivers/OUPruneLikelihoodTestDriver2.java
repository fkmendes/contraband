package testdrivers;

import beast.util.TreeParser;
import beast.evolution.tree.Node;
import contraband.prunelikelihood.OUPruneUtils;
import contraband.utils.GeneralUtils;
import org.apache.commons.math3.linear.*;

public class OUPruneLikelihoodTestDriver2 {
    public static void main(String[] args) {
        int nTraits = 4;

        double[][] h = new double[4][4];
        h[0][0] = 10;
        h[0][1] = 11;
        h[0][2] = 12;
        h[0][3] = 13;
        h[1][0] = 14;
        h[1][1] = 15;
        h[1][2] = 16;
        h[1][3] = 17;
        h[2][0] = 18;
        h[2][1] = 19;
        h[2][2] = 20;
        h[2][3] = 21;
        h[3][0] = 22;
        h[3][1] = 23;
        h[3][2] = 24;
        h[3][3] = 25;
        RealMatrix H = new Array2DRowRealMatrix(h);
        EigenDecomposition decompH = new EigenDecomposition(H);
        RealMatrix pMat = decompH.getV();
        RealMatrix inverseP = new LUDecomposition(pMat).getSolver().getInverse();


        String treeStr = "((A:23.0058179,B:23.0058179):14.350951,C:37.3567689);";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        Node aNode = tree.getNode(1);

        double[][] sigma2DArray = {{9.404673, -9.088870, 7.848381, 9.1366669},
                {0.000000, 5.281055, 1.028700, -0.9333169},
                {0.000000, 0.000000, 4.566147, 3.5514127},
                {0.000000, 0.000000, 0.000000, 5.7263340}};
        RealMatrix sigmaRM = new Array2DRowRealMatrix(sigma2DArray);

        double[][] sigmaE2DArray = {{9.404673, -9.088870, 7.848381, 9.1366669},
                {0.000000, 5.281055, 1.028700, -0.9333169},
                {0.000000, 0.000000, 4.566147, 3.5514127},
                {0.000000, 0.000000, 0.000000, 5.7263340}};
        RealMatrix sigmaERM = new Array2DRowRealMatrix(sigmaE2DArray);


        RealMatrix ouVarianceRM = OUPruneUtils.getOUVarianceRM(aNode, sigmaRM, sigmaERM, pMat, inverseP, decompH, nTraits);
        GeneralUtils.displayRealMatrix(ouVarianceRM);
    }
}
