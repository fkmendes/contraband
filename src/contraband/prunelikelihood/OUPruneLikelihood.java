package contraband.prunelikelihood;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.*;
import beast.evolution.tree.Tree;
import contraband.utils.PruneLikelihoodUtils;
import contraband.valuewrappers.OneValueContTraits;
import org.apache.commons.math3.linear.*;
import outercore.parameter.KeyRealParameter;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

@Description("This class implements likelihood for continuous traits under Ornstein–Uhlenbeck process.\n" +
        "The calculation uses Venelin's PCM likelihood.")


public class OUPruneLikelihood extends Distribution {
    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmaValuesInput = new Input<>("sigma","Trait rate matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmaeValuesInput = new Input<>("sigmae","Population error.", Input.Validate.OPTIONAL);
    final public Input<RealParameter> alphaInput = new Input<>("alpha","An array of (nTraits * (nTraits - 1) / 2) elements, representing selection strength, off-diagonal elements in Alpha matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> thetaInput = new Input<>("theta","An array of nTraits elements, representing optimum trait values, elements in Theta vector.", Input.Validate.REQUIRED);
    final public Input<KeyRealParameter> traitsValuesInput = new Input<>("traits","Trait values at tips.", Input.Validate.REQUIRED);
    final public Input<RealParameter> rootValuesInput = new Input<>("root","Trait values at the root.");

    private RealMatrix sigmaRM;
    private RealMatrix sigmaERM;
    private RealMatrix alphaRM;
    private RealVector thetaVec;
    private RealVector rootValuesVec;
    private RealMatrix identity;

    private double [] sigma;
    private double [] sigmae;
    private double[] alpha;
    private double[] theta;

    private Tree tree;
    private int nTraits;

    private double[] traitValuesArr;
    private List<RealMatrix> lMatList;
    private List<RealVector>  mVecList;
    private double [] rArr;
    private double [] vcvMatDetArr;
    private double [] negativeTwoAplusLDetArr;


    @Override
    public void initAndValidate() {
        //super.initAndValidate();
        tree = treeInput.get();
        int nodeCount = tree.getNodeCount();

        KeyRealParameter traitsValues = traitsValuesInput.get();

        nTraits = traitsValues.getMinorDimension1();
        int nSpecies = traitsValues.getMinorDimension2();

        traitValuesArr = new double[nSpecies * nTraits];
        PruneLikelihoodUtils.populateTraitValuesArr(traitsValues, tree, nTraits, traitValuesArr);

        sigmaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        populateSigmaMatrix(sigmaRM, sigmaValuesInput.get().getDoubleValues());

        alphaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        populateAlphaMatrix(alphaRM, alphaInput.get().getDoubleValues());

        thetaVec = new ArrayRealVector(new double[nTraits]);
        populateRealVector(thetaVec, thetaInput.get().getDoubleValues());

        if (sigmaeValuesInput.get() != null) {
            sigmaERM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
            populateSigmaMatrix(sigmaERM, sigmaeValuesInput.get().getDoubleValues());
        }

        lMatList = new ArrayList<>(nodeCount);
        mVecList = new ArrayList<>(nodeCount);
        rArr = new double[nodeCount];
        rootValuesVec = new ArrayRealVector(rootValuesInput.get().getDoubleValues());
        vcvMatDetArr = new double[nodeCount];
        negativeTwoAplusLDetArr = new double[nodeCount];

        // initiate
        RealMatrix iniLMat = new Array2DRowRealMatrix(new double [nTraits][nTraits]);
        RealVector iniMVec = new ArrayRealVector(new double [nTraits]);
        for (int i= 0; i < tree.getNodeCount(); i ++) {
            lMatList.add(i, iniLMat);
            mVecList.add(i, iniMVec);
            rArr[i] = 0;

        }

        identity = MatrixUtils.createRealIdentityMatrix(nTraits);
    }

    @Override
    public double calculateLogP() {


        EigenDecomposition decompositionH = new EigenDecomposition(alphaRM);
        RealMatrix pMat = decompositionH.getV();
        RealMatrix inverseP = new LUDecomposition(pMat).getSolver().getInverse();

        OUPruneUtils.pruneOUPCM(tree.getRoot(), nTraits, traitValuesArr,
                lMatList,  mVecList,  rArr,
                sigmaRM, sigmaERM,
                thetaVec, alphaRM,
                pMat, inverseP, decompositionH,
                identity, vcvMatDetArr,  negativeTwoAplusLDetArr);

        RealMatrix l0Mat = lMatList.get(tree.getRoot().getNr());
        RealVector m0Vec = mVecList.get(tree.getRoot().getNr());
        double r0 = rArr[tree.getRoot().getNr()];

        // calculate likelihood
        // loglik <- X0 * L0 * X0+ m0 * X0 + r0
        logP = l0Mat.preMultiply(rootValuesVec).dotProduct(rootValuesVec) + rootValuesVec.dotProduct(m0Vec) + r0;

        return logP;
    }

    public void populateSigmaMatrix(RealMatrix rm, double[] values) {
        int k = 0;
        for (int i = 0; i < rm.getColumnDimension(); i++) {
            for (int j = i; j < rm.getColumnDimension(); j++) {
                rm.setEntry(i, j, values[k]);
                k++;
            }
        }
    }

    public void populateAlphaMatrix(RealMatrix rm, double[] values) {
        int k = 0;
        for (int i = 0; i < rm.getColumnDimension(); i++) {
            for (int j = 0; j < rm.getColumnDimension(); j++) {
                rm.setEntry(i, j, values[k]);
                k++;
            }
        }
    }

    public void populateRealVector(RealVector vec, double[] values) {
        for (int i = 0; i < values.length; i ++) {
            vec.setEntry(i, values[i]);
        }
    }

    /*
     * This method converts trait values to
     * a list of real vectors, each real vector
     * contains all traits values of a species.
     */
    public static void populateTraitsValuesList (OneValueContTraits traitValues, Tree tree, List<RealVector> traitsValuesList) {
        int nTips = traitValues.getNSpp();
        for (int i = 0; i < nTips; i ++) {
            traitsValuesList.add(i, new ArrayRealVector(traitValues.getSpValues(tree.getNode(i).getID())));
        }
    }

    @Override
    public List<String> getArguments() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public List<String> getConditions() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        // TODO Auto-generated method stub

    }

}
