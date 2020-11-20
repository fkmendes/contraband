package contraband.prunelikelihood;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import org.apache.commons.math3.linear.*;
import outercore.parameter.KeyRealParameter;
import beast.evolution.tree.Node;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class OUNodeMath extends CalculationNode {
    final public Input<KeyRealParameter> traitsValuesInput = new Input<>("traits","Trait values at tips.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmaValuesInput = new Input<>("sigma","Trait rate matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmaeValuesInput = new Input<>("sigmae","Population error.", Input.Validate.OPTIONAL);
    final public Input<RealParameter> alphaInput = new Input<>("alpha","An array of (nTraits * (nTraits - 1) / 2) elements, representing selection strength, off-diagonal elements in Alpha matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> thetaInput = new Input<>("theta","An array of nTraits elements, representing optimum trait values, elements in Theta vector.", Input.Validate.REQUIRED);
    final public Input<RealParameter> rootValuesInput = new Input<>("root","Trait values at the root.");

    private RealMatrix sigmaRM;
    private RealMatrix sigmaERM;
    private RealMatrix alphaRM;
    private RealVector thetaVec;
    private RealVector rootValuesVec;
    private RealMatrix identity;
    private int nTraits;
    private int nSpecies;
    private int nodeNr;

    private List<RealMatrix> lMatList;
    private List<RealVector>  mVecList;
    private double [] rArr;
    private double vcvMatDet;

    private RealMatrix pMat;
    private RealMatrix inverseP;
    private EigenDecomposition decompositionH;

    private RealMatrix iniMatrix;
    private RealVector iniVec;
    private boolean singularMatrix;
    private RealMatrix invVCVMat;
    private double varianceRMDet;
    private RealMatrix aPlusLInv;
    private double logDetVNode;

    @Override
    public void initAndValidate() {
        // collect trait information
        KeyRealParameter traitsValues = traitsValuesInput.get();
        nTraits = traitsValues.getMinorDimension1();
        nSpecies = traitsValues.getMinorDimension2();

        nodeNr = 2 * nSpecies - 1;

        sigmaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        OUPruneUtils.populateSigmaMatrix(sigmaRM, sigmaValuesInput.get().getDoubleValues());
        sigmaRM = sigmaRM.multiply(sigmaRM.transpose());

        alphaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        OUPruneUtils.populateAlphaMatrix(alphaRM, alphaInput.get().getDoubleValues());

        thetaVec = new ArrayRealVector(new double[nTraits]);
        OUPruneUtils.populateRealVector(thetaVec, thetaInput.get().getDoubleValues());

        if (sigmaeValuesInput.get() != null) {
            sigmaERM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
            OUPruneUtils.populateSigmaMatrix(sigmaERM, sigmaeValuesInput.get().getDoubleValues());
            sigmaERM = sigmaERM.multiply(sigmaERM.transpose());
        }

        lMatList = new ArrayList<>(nodeNr);
        mVecList = new ArrayList<>(nodeNr);
        rArr = new double[nodeNr];
        rootValuesVec = new ArrayRealVector(rootValuesInput.get().getDoubleValues());

        // initiate
        iniMatrix = new Array2DRowRealMatrix(new double [nTraits][nTraits]);
        iniVec = new ArrayRealVector(new double [nTraits]);
        for (int i= 0; i < nodeNr; i ++) {
            lMatList.add(i, iniMatrix);
            mVecList.add(i, iniVec);
            rArr[i] = 0;

        }

        identity = MatrixUtils.createRealIdentityMatrix(nTraits);
    }

    // getters
    public RealMatrix getIniMatrix () {return iniMatrix;}

    public RealVector getIniVec () {return iniVec;}

    public RealMatrix getInverseVarianceMatrix () {return invVCVMat;}

    public double getVCVMatDet () {return  vcvMatDet;}

    public RealMatrix getLMatForNode (int nodeIdx) {return lMatList.get(nodeIdx); }

    public RealVector getMVecForNode(int nodeIdx) {return mVecList.get(nodeIdx); }

    public double getRForNode (int nodeIdx) { return rArr[nodeIdx]; }

    public RealVector getRootValuesVec () { return rootValuesVec; }

    public RealMatrix getAlphaMatrix () {return alphaRM; }

    public RealVector getThetaVec () {return thetaVec; }

    public RealMatrix getIdentityMatrix ( ) {return identity; }

    public RealMatrix getAPlusLInv() { return aPlusLInv; }

    public double getNegativeTwoAPlusLDet () { return logDetVNode; }

    public boolean getSingularMatrix() { return singularMatrix; }

    // setters
    public void setLMatForNode (int nodeIdx, RealMatrix value) { lMatList.set(nodeIdx,value); }

    public void setMVecForNode (int nodeIdx, RealVector value) { mVecList.set(nodeIdx, value); }

    public void setRForNode (int nodeIdx, double value) { rArr[nodeIdx] = value; }

    public void setSingularMatrix(boolean value) { singularMatrix = value; }


    public void populateAlphaMatrix() {
        OUPruneUtils.populateAlphaMatrix(alphaRM, alphaInput.get().getDoubleValues());
    }

    public void performAlphaDecompostion () {
        decompositionH = new EigenDecomposition(alphaRM);
        pMat = decompositionH.getV();
        inverseP = new LUDecomposition(pMat).getSolver().getInverse();
    }

    public void populateVarianceCovarianceMatrix(Node node){
        double branchLength = node.getLength();
        int nodeIdx = node.getNr();

        // variance-covariance matrix
        RealMatrix variance = calculateVarianceForNode(branchLength);

        // considering population variances
        if (node.isLeaf() && sigmaERM!=null) {
            variance = variance.add(sigmaERM);
        }

        // calculate the inverse of variance-covariance matrix
        // and determinant of variance covariance matrix
        populateInverseVarianceCovarianceMatrix(variance);
    }

    public void performAPlusOperations (Node aNode, RealMatrix aMat) {
        // (aMat + lMat).inverse
        RealMatrix AplusL = aMat.add(lMatList.get(aNode.getNr()));

        // ensure symmetry: AplusL <- 0.5 * (AplusL + AplusL.transpose)
        AplusL = (AplusL.add(AplusL.transpose())).scalarMultiply(0.5);

        try {
            // log(det(-2*AplusL))
            // determinant of -2 * (aMat + lMat) in log space
            logDetVNode = Math.log(new LUDecomposition(AplusL.scalarMultiply(-2)).getDeterminant());

            // inverse of AplusL
            // AplusL <- AplusL.inverse()
            aPlusLInv = new LUDecomposition(AplusL).getSolver().getInverse();
        }
        catch (SingularMatrixException e) {
            singularMatrix = true;
        }

        if (logDetVNode == 0.0) {
            singularMatrix = true;
        }
    }

    private RealMatrix calculateVarianceForNode(double nodeBranchLength) {
        // P_1SigmaP_t = inverseP * Sigma * t(inverseP)
        RealMatrix P_1SigmaP_t = inverseP.multiply(sigmaRM).multiply(inverseP.transpose());

        // fLambda_ij(t) * P_1SigmaP_t
        for (int i = 0; i < nTraits; i++) {
            for (int j = 0; j < nTraits; j++) {
                if (decompositionH.getRealEigenvalue(i) == 0.0 && decompositionH.getRealEigenvalue(j) == 0.0) {
                    //fLambdaP_1SigmaP_t[i][j] = nodeHeight * P_1SigmaP_t.getEntry(i,j);
                    P_1SigmaP_t.setEntry(i, j, nodeBranchLength * P_1SigmaP_t.getEntry(i, j));
                } else {
                    double lambda = decompositionH.getRealEigenvalue(i) + decompositionH.getRealEigenvalue(j);
                    double f = (1 - Math.exp(-lambda * nodeBranchLength)) / lambda;
                    //fLambdaP_1SigmaP_t[i][j] = f * P_1SigmaP_t.getEntry(i,j);
                    P_1SigmaP_t.setEntry(i, j, f * P_1SigmaP_t.getEntry(i, j));
                }
            }
        }

        // variance matrix = P * (fLambda * P_1SigmaP_t) * t(P)
        return pMat.multiply(P_1SigmaP_t).multiply(pMat.transpose());
    }

    private void populateInverseVarianceCovarianceMatrix(RealMatrix variance){

        // evaluate if nearly singular
        double[] singularValues = new SingularValueDecomposition(variance).getSingularValues();
        double min = Arrays.stream(singularValues).min().getAsDouble();
        double max = Arrays.stream(singularValues).max().getAsDouble();

        EigenDecomposition decomposition = new EigenDecomposition(variance);
        double[] eValues = decomposition.getRealEigenvalues();

        for (double ei : eValues) {
            if (ei < 1.0E-5) {
                singularMatrix = true;
            }
        }

        if ((min / max) < 1.0E-6) {
            singularMatrix = true;
        }

        // inverse of V and determinant of V
        // V <- V.inverse()
        try {
            LUDecomposition VMatLUD = new LUDecomposition(variance);
            invVCVMat = VMatLUD.getSolver().getInverse();
            vcvMatDet = VMatLUD.getDeterminant();
        } catch (SingularMatrixException e) {
            singularMatrix = true;
        }

        if (vcvMatDet == 0.0) {
            singularMatrix = true;
        }
    }

}
