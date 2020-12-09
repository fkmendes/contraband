package contraband.prunelikelihood;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
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
    
    final public Input<RealParameter> alphaInput = new Input<>("alpha","An array of (nTraits * (nTraits - 1) / 2) elements, representing selection strength, off-diagonal elements in Alpha matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> thetaInput = new Input<>("theta","An array of nTraits elements, representing optimum trait values, elements in Theta vector.", Input.Validate.REQUIRED);
    final public Input<RealParameter> rootValuesInput = new Input<>("root","Trait values at the root.");

    final public Input<RealParameter> sigmasqInput = new Input<>("sigmasq", "Evolutionary rates of traits. Diagonal elements in rate matrix.");
    final public Input<RealParameter> rhoInput = new Input<>("correlation", "Correlations of traits. Off-diagonal elements in rate matrix.");
    final public Input<RealParameter> covarianceInput = new Input<>("covariance", "cov_ij = sigma_i * sigma_j * rho_ij.");

    final public Input<RealParameter> popVarInput = new Input<>("popVar", "Variance within population. Diagonal elements in SigmaE matrix.");
    final public Input<RealParameter> popCovInput = new Input<>("popCov", "Covariance within population. Off-diagonal elements in SigmaE matrix.");

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
    private double likForSA;

    private int optNr;
    private List<RealVector> thetaVecList;
    @Override
    public void initAndValidate() {
        // collect trait information
        KeyRealParameter traitsValues = traitsValuesInput.get();
        nTraits = traitsValues.getMinorDimension1();
        nSpecies = traitsValues.getMinorDimension2();

        // number of nodes in the species tree including root
        nodeNr = 2 * nSpecies - 1;

        // if sigma matrix are specified by upper matrix
        // TRUE Sigma = sigma %*% t(sigma)
        // FALSE Sigma = 0.5 * (sigma + t(sigma))
        useUpperMatrix = useUpperMatrixInput.get();

        // Sigma matrix has trait evolutionary rates and traits correlations
        sigmaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        // check dimension
        if (sigmasqInput.get().getDimension() != nTraits) {
            sigmasqInput.get().setDimension(nTraits);
        }
        if(nTraits > 1) {
            if (rhoInput.get() == null && covarianceInput.get() != null) {
                if (covarianceInput.get().getDimension() != nTraits * (nTraits - 1) / 2) {
                    covarianceInput.get().setDimension(nTraits * (nTraits - 1) / 2);
                }
            }

            if (rhoInput.get() != null && covarianceInput.get() == null) {
                if (rhoInput.get().getDimension() != nTraits * (nTraits - 1) / 2) {
                    rhoInput.get().setDimension(nTraits * (nTraits - 1) / 2);
                }
            }
        }
        // considering population variance
        // SigmaE = sigmae %*% t(sigmae)
        if (popVarInput.get() != null) {
            sigmaERM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
            // check dimensions
            if (popVarInput.get().getDimension() != nSpecies) {
                popVarInput.get().setDimension(nTraits);
            }
            if (popCovInput.get().getDimension() != nTraits * (nTraits - 1) / 2) {
                popCovInput.get().setDimension(nTraits * (nTraits - 1) / 2);
            }
        }

        // alpha matrix is specified by nTraits * nTraits elements
        alphaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        if (alphaInput.get().getDimension() != nTraits * nTraits){
            alphaInput.get().setDimension(nTraits * nTraits);
        }

        /*
         *  One global theta vector for the tree
         */
        //thetaVec = new ArrayRealVector(new double[nTraits]);
        //OUPruneUtils.populateRealVector(thetaVec, i, thetaInput.get().getDoubleValues());

        // there are optNr theta vectors in the tree
        // corresponding to categories assignment of each branch
        optNr = optNrInput.get();
        IntegerParameter rateCategories = optAssignInput.get();
        if(optAssignInput.get().getDimension() != nodeNr) {
            rateCategories.setDimension(nodeNr);
        }
        // a list of theta vectors
        // including optNr different theta vectors
        // the index corresponds to the rate category
        thetaVecList = new ArrayList<>(optNr);
        thetaInput.get().setDimension(nTraits * optNr);
        for (int j = 0; j < optNr; j ++) {
            RealVector nodeTheta = new ArrayRealVector(new double[nTraits]);
            for (int i = 0; i < nTraits; i ++) {
                nodeTheta.setEntry(i, thetaInput.get().getDoubleValues()[j*nTraits + i]);
            }
            thetaVecList.add(j, nodeTheta);
        }

        // traits values at the root
        rootValuesInput.get().setDimension(nTraits);
        rootValuesVec = new ArrayRealVector(new double[nTraits]);

        // initialize lists
        iniMatrix = new Array2DRowRealMatrix(new double [nTraits][nTraits]);
        iniVec = new ArrayRealVector(new double [nTraits]);
        lMatList = new ArrayList<>(nodeNr);
        mVecList = new ArrayList<>(nodeNr);
        rArr = new double[nodeNr];
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
    public RealMatrix getVCVMatForNode(int nodeIdx) { return VCVMatList.get(nodeIdx);}
    public RealMatrix getInvVCVMatForNode (int nodeIdx) {return invVCVMatList.get(nodeIdx); }

    public RealMatrix getPMat () {return pMat;}

    public RealMatrix getInverseP () {return inverseP;}

    public EigenDecomposition getAlphaDecomposition() {return decompositionH;}

    public RealMatrix getVarianceMatrix() {return variance; }

    public double getLikelihoodForSA() {return likForSA;}

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
    public void setLikelihoodForSA (double value) { likForSA = value; }


    public void populateAlphaMatrix() {
        OUPruneUtils.populateAlphaMatrix(alphaRM, alphaInput.get().getDoubleValues());
    }

    public void performAlphaDecomposition() {
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
        // calculate the inverse matrix of pMat
        try {
            inverseP = new LUDecomposition(pMat).getSolver().getInverse();
        }
        catch (SingularMatrixException e) {
            singularMatrix = true;
        }
    }

    public void populateVarianceCovarianceMatrix(Node node){
        double branchLength = node.getLength();
        // deal with zero branch length node
        if(branchLength < 1e-6){
            branchLength = node.getParent().getLength();
        }
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

    public void updateSigmaMatrix(){
        if(nTraits == 1) {
            sigmaRM.setEntry(0, 0, sigmasqInput.get().getValue());
        } else {
            if (rhoInput.get() == null && covarianceInput.get() != null) {
                if (useUpperMatrix) {
                    OUPruneUtils.populateUpperSigmaMatrix(sigmaRM, sigmasqInput.get().getDoubleValues(), covarianceInput.get().getDoubleValues(), nTraits);
                    sigmaRM = sigmaRM.multiply(sigmaRM.transpose());
                } else {
                    OUPruneUtils.populateDirectSigmaMatrix(sigmaRM, sigmasqInput.get().getDoubleValues(), covarianceInput.get().getDoubleValues(), nTraits);
                }
            }

            if (rhoInput.get() != null && covarianceInput.get() == null) {
                if (useUpperMatrix) {
                    OUPruneUtils.populateUpperRhoSigmaMatrix(sigmaRM, sigmasqInput.get().getDoubleValues(), rhoInput.get().getDoubleValues(), nTraits);
                    sigmaRM = sigmaRM.multiply(sigmaRM.transpose());
                } else {
                    OUPruneUtils.populateSigmaMatrix(sigmaRM, sigmasqInput.get().getDoubleValues(), rhoInput.get().getDoubleValues(), nTraits);
                }
            }
            // considering population variance
            // SigmaE = sigmae %*% t(sigmae)
            if (popVarInput.get() != null) {
                OUPruneUtils.populateUpperSigmaMatrix(sigmaERM, popVarInput.get().getDoubleValues(), popCovInput.get().getDoubleValues(), nTraits);
                sigmaERM = sigmaERM.multiply(sigmaERM.transpose());
            }
        }
    }

    public void updateThetaVectors(){
        for (int j = 0; j < optNr; j ++) {
            RealVector nodeTheta = new ArrayRealVector(new double[nTraits]);
            for (int i = 0; i < nTraits; i ++) {
                nodeTheta.setEntry(i, thetaInput.get().getDoubleValues()[j*nTraits + i]);
            }
            thetaVecList.set(j, nodeTheta);
        }
    }

    public void updateRootValues(){
        OUPruneUtils.populateRealVector(rootValuesVec, nTraits, rootValuesInput.get().getDoubleValues());
    }
}
