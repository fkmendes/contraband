package contraband.prunelikelihood;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import org.apache.commons.math3.exception.MathUnsupportedOperationException;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;
import outercore.parameter.KeyRealParameter;
import beast.evolution.tree.Node;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class OUNodeMath extends CalculationNode {
    final public Input<KeyRealParameter> traitsValuesInput = new Input<>("traits","Trait values at tips.", Input.Validate.REQUIRED);
    final public Input<Boolean> useUpperMatrixInput = new Input<>("upperMatrix", "TRUE, if sigmasq and correlations are from upper matrix", true);
    final public Input<RealParameter> sigmaValuesInput = new Input<>("sigma","Trait rate matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmaeValuesInput = new Input<>("sigmae","Population error.", Input.Validate.OPTIONAL);
    final public Input<RealParameter> alphaInput = new Input<>("alpha","An array of (nTraits * (nTraits - 1) / 2) elements, representing selection strength, off-diagonal elements in Alpha matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> thetaInput = new Input<>("theta","An array of nTraits elements, representing optimum trait values, elements in Theta vector.", Input.Validate.REQUIRED);
    final public Input<RealParameter> rootValuesInput = new Input<>("root","Trait values at the root.");

    final public Input<Integer> optNrInput = new Input<>("optNr","Number of theta (vectors).");
    final public Input<IntegerParameter> optAssignInput = new Input<>("optAssign", "the opt assignment for each node in the tree.");

    private Boolean useUpperMatrix;
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

    private List<RealMatrix> AMatList;
    private List<RealMatrix> CMatList;
    private List<RealMatrix> EMatList;
    private List<RealVector>  bVecList;
    private List<RealVector>  dVecList;
    private double [] fArr;

    private List<RealMatrix> PhiMatList;
    private List<RealVector>  omegaList;

    private List<RealMatrix> VCVMatList;
    private List<RealMatrix> invVCVMatList;

    private double vcvMatDet;

    private RealMatrix pMat;
    private RealMatrix inverseP;
    private EigenDecomposition decompositionH;

    private RealMatrix iniMatrix;
    private RealMatrix variance;
    private RealVector iniVec;
    private boolean singularMatrix;
    private RealMatrix invVCVMat;
    private RealMatrix aPlusLInv;
    private double logDetVNode;

    private int optNr;
    private List<RealVector> thetaVecList;
    @Override
    public void initAndValidate() {
        // collect trait information
        KeyRealParameter traitsValues = traitsValuesInput.get();
        nTraits = traitsValues.getMinorDimension1();
        nSpecies = traitsValues.getMinorDimension2();

        nodeNr = 2 * nSpecies - 1;

        useUpperMatrix = useUpperMatrixInput.get();

        sigmaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        if (useUpperMatrix) {
            OUPruneUtils.populateSigmaMatrix(sigmaRM, sigmaValuesInput.get().getDoubleValues());
            sigmaRM = sigmaRM.multiply(sigmaRM.transpose());
        } else {
            double[] sigmas = sigmaValuesInput.get().getDoubleValues();
            int k = 0;
            for (int i = 0; i < nTraits; i++) {
                for (int j = i; j < nTraits; j++) {
                    sigmaRM.setEntry(i, j, sigmas[k]);
                    sigmaRM.setEntry(j, i, sigmas[k]);
                    k++;
                }
            }
        }

        alphaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        OUPruneUtils.populateAlphaMatrix(alphaRM, alphaInput.get().getDoubleValues());

        thetaVec = new ArrayRealVector(new double[nTraits]);
        optNr = optNrInput.get();
        thetaVecList = new ArrayList<>(optNr);
        for (int j = 0; j < optNr; j ++) {
            //OUPruneUtils.populateRealVector(thetaVec, i, thetaInput.get().getDoubleValues());
            RealVector nodeTheta = new ArrayRealVector(new double[nTraits]);
            for (int i = 0; i < nTraits; i ++) {
                nodeTheta.setEntry(i, thetaInput.get().getDoubleValues()[j*nTraits + i]);
            }
            thetaVecList.add(j, nodeTheta);
        }

        if (sigmaeValuesInput.get() != null) {
            sigmaERM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
            OUPruneUtils.populateSigmaMatrix(sigmaERM, sigmaeValuesInput.get().getDoubleValues());
            sigmaERM = sigmaERM.multiply(sigmaERM.transpose());
        }

        lMatList = new ArrayList<>(nodeNr);
        mVecList = new ArrayList<>(nodeNr);
        rArr = new double[nodeNr];
        rootValuesVec = new ArrayRealVector(rootValuesInput.get().getDoubleValues());
        fArr = new double[nodeNr];
        AMatList = new ArrayList<>(nodeNr);
        CMatList = new ArrayList<>(nodeNr);
        EMatList = new ArrayList<>(nodeNr);
        PhiMatList = new ArrayList<>(nodeNr);
        VCVMatList = new ArrayList<>(nodeNr);
        invVCVMatList = new ArrayList<>(nodeNr);
        bVecList = new ArrayList<>(nodeNr);
        dVecList = new ArrayList<>(nodeNr);
        omegaList = new ArrayList<>(nodeNr);
        // initiate
        iniMatrix = new Array2DRowRealMatrix(new double [nTraits][nTraits]);
        iniVec = new ArrayRealVector(new double [nTraits]);
        for (int i= 0; i < nodeNr; i ++) {
            lMatList.add(i, iniMatrix);
            mVecList.add(i, iniVec);
            AMatList.add(i, iniMatrix);
            CMatList.add(i, iniMatrix);
            EMatList.add(i, iniMatrix);
            bVecList.add(i, iniVec);
            dVecList.add(i, iniVec);
            PhiMatList.add(i, iniMatrix);
            omegaList.add(i, iniVec);
            VCVMatList.add(i, iniMatrix);
            invVCVMatList.add(i, iniMatrix);
        }

        identity = MatrixUtils.createRealIdentityMatrix(nTraits);
        variance = new Array2DRowRealMatrix(new double [nTraits][nTraits]);
        pMat= new Array2DRowRealMatrix(new double [nTraits][nTraits]);

        IntegerParameter rateCategories = optAssignInput.get();
        if(optAssignInput.get().getDimension() != nodeNr) {
            rateCategories.setDimension(nodeNr);
        }

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

    public RealVector getThetaForNode(int nodeIdx) {
        return thetaVecList.get(optAssignInput.get().getValues()[nodeIdx]);
    }

    public RealMatrix getAMatForNode (int nodeIdx) {return AMatList.get(nodeIdx);}
    public RealMatrix getCMatForNode (int nodeIdx) {return CMatList.get(nodeIdx);}
    public RealMatrix getEMatForNode (int nodeIdx) {return EMatList.get(nodeIdx);}
    public RealMatrix getPhiForNode (int nodeIdx) {return PhiMatList.get(nodeIdx);}
    public RealVector getbVecForNode (int nodeIdx) {return bVecList.get(nodeIdx);}
    public RealVector getdVecForNode (int nodeIdx) {return dVecList.get(nodeIdx);}
    public RealVector getOmegaForNode (int nodeIdx) {return omegaList.get(nodeIdx);}
    public double getfForNode(int nodeIdx) {return fArr[nodeIdx];}

    public RealMatrix getPMat () {return pMat;}

    public RealMatrix getInverseP () {return inverseP;}

    public EigenDecomposition getAlphaDecomposition() {return decompositionH;}

    public RealMatrix getVarianceMatrix() {return variance; }

    // setters
    public void setLMatForNode (int nodeIdx, RealMatrix value) { lMatList.set(nodeIdx,value); }

    public void setMVecForNode (int nodeIdx, RealVector value) { mVecList.set(nodeIdx, value); }

    public void setRForNode (int nodeIdx, double value) { rArr[nodeIdx] = value; }

    public void setSingularMatrix(boolean value) { singularMatrix = value; }

    public void setAbCdEfOmegaPhiForNode (int nodeIdx, RealMatrix A, RealVector b, RealMatrix C, RealVector d, RealMatrix E, double f, RealVector omega, RealMatrix Phi) {
        AMatList.set(nodeIdx,A);
        CMatList.set(nodeIdx,C);
        EMatList.set(nodeIdx,E);
        PhiMatList.set(nodeIdx,Phi);
        bVecList.set(nodeIdx, b);
        dVecList.set(nodeIdx, d);
        omegaList.set(nodeIdx, omega);
        fArr[nodeIdx] = f;

    }


    public void populateAlphaMatrix() {
        OUPruneUtils.populateAlphaMatrix(alphaRM, alphaInput.get().getDoubleValues());
    }

    public void performAlphaDecompostion () {
        decompositionH = new EigenDecomposition(alphaRM);
        // normalize eigen vectors
        // NOTE: in R, eigen() returns vectors in decreasing order according to eigen values
        for (int i = 0; i < nTraits; i ++) {
            RealVector v = decompositionH.getEigenvector(i);
            double sum = 0.0;
            for (int j = 0; j < nTraits; j++){
                sum += (v.getEntry(j) * v.getEntry(j));
            }
            pMat.setColumnVector(i,v.mapMultiply(1.0/FastMath.sqrt(sum)));
        }

        inverseP = new LUDecomposition(pMat).getSolver().getInverse();
    }

    public void populateVarianceCovarianceMatrix(Node node){
        double branchLength = node.getLength();
        int nodeIdx = node.getNr();

        // variance-covariance matrix
        variance = calculateVarianceForNode(branchLength);

        // considering population variances
        if (node.isLeaf() && sigmaERM!=null) {
            variance = variance.add(sigmaERM);
        }
        // set variance
        VCVMatList.set(nodeIdx, variance);

        // calculate the inverse of variance-covariance matrix
        // and determinant of variance covariance matrix
        populateInverseVarianceCovarianceMatrix(variance);
        invVCVMatList.set(nodeIdx, invVCVMat);
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
        RealMatrix fLambdaP_1SigmaP_t = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        for (int i = 0; i < nTraits; i++) {
            for (int j = 0; j < nTraits; j++) {
                double lambda = decompositionH.getRealEigenvalue(i) + decompositionH.getRealEigenvalue(j);
                if (Math.abs(lambda) < 1e-8) {
                    //fLambdaP_1SigmaP_t[i][j] = nodeHeight * P_1SigmaP_t.getEntry(i,j);
                    fLambdaP_1SigmaP_t.setEntry(j, i, P_1SigmaP_t.getEntry(i,j) * nodeBranchLength);
                } else {
                    double f = (1 - Math.exp(-lambda * nodeBranchLength)) / lambda;
                    //fLambdaP_1SigmaP_t[i][j] = f * P_1SigmaP_t.getEntry(i,j);
                    fLambdaP_1SigmaP_t.setEntry(j, i, f * P_1SigmaP_t.getEntry(i,j));
                }
            }
        }

        // variance matrix = P * (fLambda * P_1SigmaP_t) * t(P)
        return pMat.multiply(fLambdaP_1SigmaP_t).multiply(pMat.transpose());
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
