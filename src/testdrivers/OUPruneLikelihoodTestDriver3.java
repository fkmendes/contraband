package testdrivers;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.prunelikelihood.OUPruneLikelihood;
import contraband.prunelikelihood.OUPruneUtils;
import contraband.utils.GeneralUtils;
import contraband.valuewrappers.OneValueContTraits;
import org.apache.commons.math3.linear.*;
import java.util.ArrayList;
import java.util.List;
import beast.evolution.tree.Node;

public class OUPruneLikelihoodTestDriver3 {
    public static void main(String[] args) {
        String treeStr = "((A:3.0058179,B:3.0058179):4.350951,C:7.3567689);";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        int nTraits = 2;

        RealParameter data = new RealParameter(new Double[] {
                0.131394584822684, -0.65643500381027,
                -0.17399894787897, -0.19269144948362,
                -0.373259036803319, 1.27056761078824
        });
        String spNames = "A,B,C";
        OneValueContTraits traitValues = new OneValueContTraits();
        traitValues.initByName("nTraits", nTraits, "spNames", spNames, "traitValues", data);

        // OU model parameters
        RealParameter rootValues = new RealParameter(new Double[] {0.2, 1.3});

        RealParameter sigma = new RealParameter(new Double[] {1.0, 0.6, 1.0});

        RealParameter sigmae = new RealParameter(new Double[] {1.0, 0.3, 1.0});

        RealParameter alpha = new RealParameter(new Double[] {2.0, 0.0, 0.0, 2.0});

        RealParameter theta = new RealParameter(new Double[] {0.5, 0.5});

        // TEST
        OUPruneLikelihood pcm = new OUPruneLikelihood();
        int nSpecies = tree.getLeafNodeCount();
        int nodeCount = tree.getNodeCount();


        List<RealVector> traitsValuesList = new ArrayList<>(nSpecies);
        pcm.populateTraitsValuesList(traitValues, tree, traitsValuesList);

        RealMatrix sigmaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        pcm.populateSigmaMatrix(sigmaRM, sigma.getDoubleValues());

        RealMatrix alphaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        pcm.populateAlphaMatrix(alphaRM, alpha.getDoubleValues());

        RealVector thetaVec = new ArrayRealVector(new double[nTraits]);
        pcm.populateRealVector(thetaVec, theta.getDoubleValues());

        RealMatrix sigmaERM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        pcm.populateSigmaMatrix(sigmaERM, sigmae.getDoubleValues());


        List<RealMatrix> lMatList = new ArrayList<>(nodeCount);
        List<RealVector> mVecList = new ArrayList<>(nodeCount);
        double [] rArr = new double[nodeCount];
        RealVector rootValuesVec = new ArrayRealVector(rootValues.getDoubleValues());
        double []vcvMatDetArr = new double[nodeCount];
        double [] negativeTwoAplusLDetArr = new double[nodeCount];

        // initiate
        RealMatrix iniLMat = new Array2DRowRealMatrix(new double [nTraits][nTraits]);
        RealVector iniMVec = new ArrayRealVector(new double [nTraits]);
        for (int i= 0; i < tree.getNodeCount(); i ++) {
            lMatList.add(i, iniLMat);
            mVecList.add(i, iniMVec);
            rArr[i] = 0;

        }

        RealMatrix identity = MatrixUtils.createRealIdentityMatrix(nTraits);

        EigenDecomposition decompositionH = new EigenDecomposition(alphaRM);
        RealMatrix pMat = decompositionH.getV();
        RealMatrix inverseP = new LUDecomposition(pMat).getSolver().getInverse();
        //System.out.println("Display H :");
        //GeneralUtils.displayRealMatrix(alphaRM);




        // calculate likelihood
        // loglik <- X0 * L0 * X0+ m0 * X0 + r0
        OUPruneUtils.pruneOUPCM(tree.getRoot(), nTraits, traitsValuesList, lMatList,  mVecList,  rArr, sigmaRM,  sigmaERM, thetaVec, alphaRM, pMat, inverseP, decompositionH, identity, vcvMatDetArr,  negativeTwoAplusLDetArr);
        RealMatrix l0Mat = lMatList.get(tree.getRoot().getNr());
        RealVector m0Vec = mVecList.get(tree.getRoot().getNr());
        double r0 = rArr[tree.getRoot().getNr()];
        System.out.println("Display LE:");
        GeneralUtils.displayRealMatrix(l0Mat);

        System.out.println("Display mE:");
        GeneralUtils.displayRealVector(m0Vec);

        System.out.println("Display rE = " + r0);
        double logP = l0Mat.preMultiply(rootValuesVec).dotProduct(rootValuesVec) + rootValuesVec.dotProduct(m0Vec) + r0;
        System.out.println("Log likelihood = " + logP);

        // calculate likelihood
        //OUPruneLikelihood pcm = new OUPruneLikelihood();
        //pcm.initByName("tree", tree, "traits", traitValues, "alpha", alpha, "theta", theta, "sigma", sigma, "sigmae", sigmae, "root", rootValues);
        //double loglik = pcm.calculateLogP();
        //System.out.println("Log likelihood = " + loglik);
    }
}
