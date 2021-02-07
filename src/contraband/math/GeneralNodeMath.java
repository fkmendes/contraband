package contraband.math;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import contraband.prunelikelihood.MorphologicalData;
import contraband.prunelikelihood.SigmaMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.FastMath;
import java.util.Arrays;

public class GeneralNodeMath extends CalculationNode {
    final public Input<MorphologicalData> traitInput = new Input<>("trait","Morphological data set.");
    final public Input<SigmaMatrix> rateMatrixInput = new Input<>("rateMatrix", "Trait evolutionary variance-covariance matrix, i.e. trait rate matrix.", Input.Validate.REQUIRED);
    final public Input<SigmaMatrix> popMatrixInput = new Input<>("popMatrix", "Intraspecific variance-covariance matrix, i.e. population noise.");
    final public Input<Tree> treeInput = new Input<>("tree", "Phylogenetic tree.");
    final public Input<Boolean> shareRhoInput = new Input<>("shareCorrelation", "TRUE, if correlations are shared between rateMatrix and popMatrix.", false);
    final public Input<RealParameter> rootValuesInput = new Input<>("rootValues", "Trait values at the root.");


    private Tree tree;
    private int nodeNr;

    private boolean shareRho;
    private boolean matrixParams;
    private boolean popVariance;
    private boolean sampleRoot;

    private SigmaMatrix rateMatrix;
    private SigmaMatrix popMatrix;

    private int nTraits;
    private int traitDim;
    private int paramDim;

    private double [] aArray;
    private double [] cArray;
    private double [] eArray;
    private double [] fArray;
    private double [] lArray;
    private double [] mVecArray;
    private double [] rArray;

    // temporary variables
    private double [] mVecInit;
    private double [] mVec;
    private double [] lMatInit;
    private double [] lMat;
    private double [] aMat;
    private double [] cMat;
    private double [] eMat;

    // for doing LUDecompostion
    private double[] lu;
    private int[] pivot;
    private boolean[] evenSingular;
    private double [] identityMatrix;
    private boolean singularMatrix;

    // variance and expectation for BM
    private double [] nodeVariance;
    private double [] nodeExpectation;
    private double [] expect1;
    private double [] expect2;
    private double [] expectp;

    private double[] rootValuesVec;
    private double[] storedRootValuesVec;

    // temporary variables
    private double[] traitsVec;
    private RealMatrix realMatrix;

    // this is for trees with sampled ancestors
    private double likForSA;
    private double[] traitVecForSA;

    private double [] traitRateMatrix;
    private double detTraitRateMat;
    private double [] invTraitRateMatrix;
    private double detInvTraitRateMat;

    // for matrixParams
    private double [] varianceMatrix;
    private double detVarianceMat;
    private double [] invVarianceMatrix;
    private double detInvVarianceMat;
    private double[] rateScaler;
    private double[] aPlusLInverse;
    private double logNegativeTwoAplusL;
    private double[] popVarianceMatrix;

    private double storedDetTraitRateMat;
    private double [] storedInvTraitRateMatrix;
    private double storedDetInvTraitRateMat;
    private double [] storedTraitRateMatrix;
    private double [] storedVarianceMatrix;
    private double [] storedInvVarianceMatrix;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        nodeNr = tree.getNodeCount();

        nTraits = traitInput.get().getTotalTraitNr();
        traitDim = nTraits * nTraits; // the dimension of trait rate matrix
        paramDim = 1; // the dimension of A, C, E and L , = 1 if NOT using matrix parameterization

        // if trait rate matrix and intraspecific matrix share correlation
        shareRho = shareRhoInput.get();

        // get trait rate matrix
        rateMatrix = rateMatrixInput.get();
        rateMatrix.populateSigmaMatrix();
        realMatrix = new Array2DRowRealMatrix(new double[nTraits][nTraits]);

        // get intraspecific variance-covariance matrix
        if(popMatrixInput.get() == null){
            // ignore intraspecific variance
            popVariance = false;
            matrixParams = false;
        } else {
            // account for intraspecific variance
            popVariance = true;
            popMatrix = popMatrixInput.get();
            matrixParams = !rateMatrix.getOneRateOnly() || !popMatrix.getOneRateOnly();
            if(!shareRho || matrixParams){
                if(traitInput.get().getTransformDataFlag()) {
                    throw new RuntimeException("GeneralNodeMath:: Data cannot be transformed!.");
                }
            }
            if(matrixParams) {
                // if using matrix parameterization, A, C, E and L are nTraits * nTraits matrix
                paramDim = nTraits * nTraits;
                initMatrixParams();
                popMatrix.populateSigmaMatrix();
            } else {
                popMatrix.populateSigmaValue();
            }
            popVarianceMatrix = new double[paramDim];
        }
        if(popVariance && !matrixParams) {
            // taking intraspecific variance into account
            // all traits share (1) the same evolutionary rate (2) the same trait correlations (3) the same intraspecific variance
            rateMatrix.populateRhoMatrix();
            rateMatrix.populateSigmaValue();
        }

        // get root values
        rootValuesVec = new double [nTraits];
        if(rootValuesInput.get() != null) {
            if (rootValuesInput.get().getDimension() != nTraits) {
                rootValuesInput.get().setDimension(nTraits);
            }
            rootValuesVec = rootValuesInput.get().getDoubleValues();
            sampleRoot = true;
        } else {
            Log.warning.println("NodeMath::Estimating root state by MLE.");
            sampleRoot = false;
        }
        storedRootValuesVec = new double [nTraits];
        traitsVec = new double [nTraits];


        // initialize matrix and vectors
        initAbCdEfArray();
        initLmrArray();
        initUDecomposition();
        initForSampledAncestors();
        initVarianceAndExpectation();
        initTraitRateMatrix();
    }

    // initialise
    private void initAbCdEfArray() {
        // f is a single double value for each node
        fArray = new double[nodeNr];

        // A C E are matrix for each node
        aMat = new double[paramDim];
        cMat = new double[paramDim];
        eMat = new double[paramDim];
        aArray = new double[paramDim * nodeNr];
        cArray = new double[paramDim * nodeNr];
        eArray = new double[paramDim * nodeNr];

    }

    private void initLmrArray() {
        // each nodes has a vector m
        mVecArray = new double[nTraits * nodeNr];
        mVec = new double[nTraits];
        mVecInit = new double[nTraits];
        // each nodes has a single double r
        rArray = new double[nodeNr];
        // depend on input, L can be either a matrix or a double value
        lMat = new double[paramDim];
        lMatInit = new double[paramDim];
        lArray = new double[paramDim * nodeNr];

    }

    private void initUDecomposition() {
        lu = new double [traitDim];
        pivot = new int[nTraits];
        evenSingular = new boolean[2];
        identityMatrix = new double[traitDim];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, nTraits);
        }
    }

    private void initVarianceAndExpectation() {
        // each node has its variance = branch rate * branch length * trait rate matrix + intraspecific variance (if necessary)
        if(matrixParams){
            nodeVariance = new double[nodeNr * traitDim];
        } else {
            nodeVariance = new double[nodeNr];
        }

        // each node has its expectation = weight average of trait values from descendants
        nodeExpectation = new double [nTraits * nodeNr];
        // initialize expectation
        expect1 = new double [nTraits];
        expect2 = new double [nTraits];
        expectp = new double [nTraits];
    }

    private void initForSampledAncestors () {
        likForSA = 0.0;
        traitVecForSA = new double [nTraits];
    }

    private void initTraitRateMatrix () {
        traitRateMatrix = new double[traitDim];
        invTraitRateMatrix = new double[traitDim];
        storedTraitRateMatrix = new double[traitDim];
        storedInvTraitRateMatrix = new double[traitDim];
    }

    private void initMatrixParams(){
        varianceMatrix = new double [traitDim];
        invVarianceMatrix = new double [traitDim];
        storedVarianceMatrix = new double [traitDim];
        storedInvVarianceMatrix = new double [traitDim];
        rateScaler = new double[traitDim];
        aPlusLInverse = new double [traitDim];
    }

    //getters
    public double[] getAMatForNode(int nodeIdx) {
        MatrixUtilsContra.getMatrixRow(aArray, nodeIdx, paramDim, aMat);
        return aMat;
    }

    public double[] getCMatForNode(int nodeIdx) {
        MatrixUtilsContra.getMatrixRow(cArray, nodeIdx, paramDim, cMat);
        return cMat;
    }

    public double[] getEMatForNode(int nodeIdx) {
        MatrixUtilsContra.getMatrixRow(eArray, nodeIdx, paramDim, eMat);
        return eMat;
    }

    public double getfForNode (int nodeIdx) { return fArray[nodeIdx]; }

    public double[] getLMatForNode(int nodeIdx) {
        MatrixUtilsContra.getMatrixRow(lArray, nodeIdx, paramDim, lMat);
        return lMat;
    }

    public double getRForNode (int nodeIdx) { return rArray[nodeIdx]; }

    public double[] getMVecForNode (int nodeIdx) {
        MatrixUtilsContra.getMatrixRow(mVecArray, nodeIdx, nTraits, mVec);
        return mVec;
    }

    public boolean getSingularMatrixFlag () { return singularMatrix; }

    public double getTraitRateMatrixDeterminant () { return detTraitRateMat; }

    public double[] getRootValuesArr () { return rootValuesVec; }

    public double[] getTraitRateMatrixInverse () { return invTraitRateMatrix; }

    public double getLikelihoodForSampledAncestors () { return likForSA; }

    public double[] getSampledAncestorTraitsVec () { return traitVecForSA; }

    public double[] getInitMVec () { return mVecInit; }

    public double[] getInitLMat () { return lMatInit; }

    public double[] getTraitsVec () {return traitsVec; }

    public double[] getTempVec () { return mVec; }

    public double[] getTraitRateMatrix () {
        if(popVariance && !matrixParams){
            return  rateMatrix.getRhoMatrix();
        }
        return rateMatrix.getSigmaMatrix();
    }

    public double getTraitRateMatrixInverseDeterminant () { return detInvTraitRateMat; }

    public double getVarianceForNode (int nodeIdx) { return nodeVariance[nodeIdx]; }

    public boolean getMatrixParamsFlag() { return matrixParams; }

    public boolean getPopVarianceFlag() { return popVariance; }

    public double[] getPopVarianceMatrix () { return popMatrix.getSigmaMatrix(); }

    public double getVarianceMatrixDet() { return detVarianceMat; }

    public double[] getInvVarianceMatrix() {return invVarianceMatrix; }

    public double[] getAPlusLInverse () { return aPlusLInverse; }

    public double getLogDetNegativeTwoAplusL() { return logNegativeTwoAplusL; }

    public double getTraitRate() { return rateMatrix.getSigmaMatrix()[0]; }

    // setters
    public void setLikelihoodForSampledAncestors(double value) {
        likForSA = value;
    }

    public void setAForNode (int nodeIdx, double value) { aArray[nodeIdx] = value; }

    public void setAMatForNode (int nodeIdx, double[] value) {
        System.arraycopy(value, 0, aArray, nodeIdx * traitDim, value.length);
    }

    public void setCForNode (int nodeIdx, double value) { cArray[nodeIdx] = value; }

    public void setCMatForNode (int nodeIdx, double[] value) {
        System.arraycopy(value, 0, cArray, nodeIdx * traitDim, value.length);
    }

    public void setEForNode (int nodeIdx, double value) { eArray[nodeIdx] = value; }

    public void setEMatForNode (int nodeIdx, double[] value) {
        System.arraycopy(value, 0, eArray, nodeIdx * traitDim, value.length);
    }

    public void setfForNode (int nodeIdx, double value) { fArray[nodeIdx] = value; }

    public void setLForNode (int nodeIdx, double[] value) {
       System.arraycopy(value, 0, lArray, nodeIdx * paramDim, paramDim);
    }

    public void setRForNode (int nodeIdx, double value) { rArray[nodeIdx] = value; }

    public void setMVecForNode (int nodeIdx, double[] value) {
        MatrixUtilsContra.setMatrixRow(mVecArray, value, nodeIdx, nTraits);
    }

    public void setTraitsVecForSampledAncestor (double[] traitValues, int nodeIdx) {
        MatrixUtilsContra.getMatrixRow(traitValues, nodeIdx, nTraits, traitVecForSA);
    }

    public void setTraitsVecForTip (double[] traitValues, int tipIdx) {
        MatrixUtilsContra.getMatrixRow(traitValues, tipIdx, nTraits, traitsVec);
    }

    public void setExpectationForTip (int nodeIdx) {
        MatrixUtilsContra.setMatrixRow(nodeExpectation, traitsVec, nodeIdx, nTraits);
    }

    public void setVarianceForTip (int nodeIdx, double value) { nodeVariance[nodeIdx] = value; }

    public void setVarianceForParent (int parentIdx, double branchLength, int child1Idx, int child2Idx) {
        double vc1 = nodeVariance[child1Idx];
        double vc2 = nodeVariance[child2Idx];
        nodeVariance[parentIdx] = branchLength + ((vc1 * vc2)/(vc1 + vc2));
    }

    public void setExpectationForParent (int parentIdx, int child1Idx, int child2Idx) {
        double vc1 = nodeVariance[child1Idx];
        double vc2 = nodeVariance[child2Idx];
        MatrixUtilsContra.getMatrixRow(nodeExpectation, child1Idx, nTraits, expect1);
        MatrixUtilsContra.getMatrixRow(nodeExpectation, child2Idx, nTraits, expect2);
        MatrixUtilsContra.vectorMapMultiply(expect1, vc2/(vc1+vc2), expect1);
        MatrixUtilsContra.vectorMapMultiply(expect2, vc1/(vc1+vc2), expect2);
        MatrixUtilsContra.vectorAdd(expect1, expect2, expectp);
        MatrixUtilsContra.setMatrixRow(nodeExpectation, expectp, parentIdx, nTraits);
    }

    public void setSingularMatrixFlag(boolean value) {
        singularMatrix  = value;
    }

    //
    public void populateRootValuesVec(int rootIdx) {
        if(sampleRoot && rootValuesInput.isDirty()) {
            rootValuesVec = rootValuesInput.get().getDoubleValues();
        }
        if(!sampleRoot) {
            MatrixUtilsContra.getMatrixRow(nodeExpectation, rootIdx, nTraits, rootValuesVec);
        }
    }

    // invert trait rate matrix and get determinant
    // NOTE: if we sample the inverse matrix, this method can be simplified.
    public void operateOnTraitRateMatrix() {
        singularMatrix = false;

        if(popVariance && !matrixParams){
            // taking intraspecific variance into account
            // all traits share (1) the same evolutionary rate (2) the same trait correlations (3) the same intraspecific variance
            // LUDecomposition
            LUDecompositionForArray.ArrayLUDecomposition(rateMatrix.getRhoMatrix(), lu, pivot, evenSingular, nTraits);
        } else {
            // LUDecomposition
            LUDecompositionForArray.ArrayLUDecomposition(rateMatrix.getSigmaMatrix(), lu, pivot, evenSingular, nTraits);
        }

        // invert the traitRateMatrix
        try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evenSingular[1], nTraits, invTraitRateMatrix);
        } catch (RuntimeException e) {
            singularMatrix = true;
        }
        // get the determinant of traitRateMatrix
        double det = LUDecompositionForArray.getDeterminant(lu, nTraits, evenSingular);
        if (det == 0.0) {
            singularMatrix = true;
        }
        detTraitRateMat = FastMath.log(det);
    }

    public void operateOnInvTraitRateMatrix() {
        // LUDecomposition
        LUDecompositionForArray.ArrayLUDecomposition(invTraitRateMatrix, lu, pivot, evenSingular, nTraits);

        // get the determinant of invTraitRateMatrix
        if(traitInput.get().getTransformDataFlag()) {
            detInvTraitRateMat = FastMath.log(LUDecompositionForArray.getDeterminant(lu, nTraits, evenSingular));
        } else {
            detInvTraitRateMat = LUDecompositionForArray.getDeterminant(lu, nTraits, evenSingular);
        }
    }

    public boolean updateSigmaMatrix(){
        boolean update1 = false;
        if(rateMatrixInput.isDirty()) {
            update1 = rateMatrix.updateParameters(!shareRho && matrixParams);
        }

        boolean update2 = false;
        if(popMatrixInput.isDirty()){
            update2 = popMatrix.updateParameters(!shareRho && matrixParams);
        }
        return update1 || update2;
    }

    /*
     * matrixParams operations
     */
    public void checkNearlySingularMatrix () {
        singularMatrix = false;
        // NOTE: if matrix parameterization
        // we need to check the variance matrix
        if(matrixParams) {
            for (int i = 0; i < nTraits; i++) {
                for (int j = 0; j < nTraits; j++) {
                    realMatrix.setEntry(i, j, varianceMatrix[i * nTraits + j]);
                }
            }
        } else {
            // otherwise we check the trait rate matrix
            for (int i = 0; i < nTraits; i++) {
                for (int j = 0; j < nTraits; j++) {
                    realMatrix.setEntry(i, j, rateMatrixInput.get().getSigmaMatrix()[i * nTraits + j]);
                }
            }
        }
        double[] singularValues = new SingularValueDecomposition(realMatrix).getSingularValues();
        double min = Arrays.stream(singularValues).min().getAsDouble();
        double max = Arrays.stream(singularValues).max().getAsDouble();

        EigenDecomposition decomposition = new EigenDecomposition(realMatrix);
        double[] eValues = decomposition.getRealEigenvalues();

        for (double ei : eValues) {
            if (ei < 1.0E-5) {
                singularMatrix = true;
            }
        }

        if ((min / max) < 1.0E-6) {
            singularMatrix = true;
        }
    }

    public void populateVarianceMatrixForTips(double branchLength){
        MatrixUtilsContra.vectorMapMultiply(rateMatrix.getSigmaMatrix(), branchLength, rateScaler);
        MatrixUtilsContra.vectorAdd(rateScaler, popMatrix.getSigmaMatrix(), varianceMatrix);
    }

    public void populateVarianceMatrixForIntNodes(double branchLength){
        MatrixUtilsContra.vectorMapMultiply(rateMatrix.getSigmaMatrix(), branchLength, varianceMatrix);
    }

    public void operateOnVarianceMatrix() {
        checkNearlySingularMatrix();

        // LUDecomposition
        LUDecompositionForArray.ArrayLUDecomposition(varianceMatrix, lu, pivot, evenSingular, nTraits);

        // invert the traitRateMatrix
        try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evenSingular[1], nTraits, invVarianceMatrix);
        } catch (RuntimeException e) {
            singularMatrix = true;
        }
        // get the determinant of traitRateMatrix
        double det = LUDecompositionForArray.getDeterminant(lu, nTraits, evenSingular);
        if (det == 0.0) {
            singularMatrix = true;
        }

        detVarianceMat = FastMath.log(det);
    }

    public void operateOnAPlusLMatrixForNode (int nodeIdx){
        // (A + L)
        double[] aPlusL = new double[nTraits * nTraits];
        MatrixUtilsContra.vectorAdd(getAMatForNode(nodeIdx), getLMatForNode(nodeIdx), aPlusL);

        // invert (A + L)
        LUDecompositionForArray.ArrayLUDecomposition(aPlusL, lu, pivot, evenSingular, nTraits);
        try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evenSingular[1], nTraits, aPlusLInverse);
        } catch (RuntimeException e) {
            singularMatrix = true;
        }

        // get the determinant of (A + L)
        double det = LUDecompositionForArray.getDeterminant(lu, nTraits, evenSingular);

        // get the determinant of -2 * (A + L) = (-2)^nTraits * det(A + L)
        double detVNode = FastMath.pow(-2, nTraits) * det;
        logNegativeTwoAplusL = FastMath.log(detVNode);
    }


    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    public void store() {
        super.store();

        storedDetTraitRateMat = detTraitRateMat;
        storedDetInvTraitRateMat = detInvTraitRateMat;

        System.arraycopy(rootValuesVec, 0, storedRootValuesVec, 0, nTraits);
        System.arraycopy(invTraitRateMatrix, 0, storedInvTraitRateMatrix, 0, traitDim);

        if(matrixParams) {
            System.arraycopy(varianceMatrix, 0, storedVarianceMatrix, 0, traitDim);
            System.arraycopy(invVarianceMatrix, 0, storedInvVarianceMatrix, 0, traitDim);
        }
    }

    @Override
    public void restore() {
        super.restore();

        double tempDetTraitRateMat = detTraitRateMat;
        detTraitRateMat = storedDetTraitRateMat;
        storedDetTraitRateMat = tempDetTraitRateMat;

        double tempDetInvTraitRateMat = detInvTraitRateMat;
        detInvTraitRateMat = storedDetInvTraitRateMat;
        storedDetInvTraitRateMat = tempDetInvTraitRateMat;

        double[] tempInvTraitRateMatrix = invTraitRateMatrix;
        invTraitRateMatrix = storedInvTraitRateMatrix;
        storedInvTraitRateMatrix = tempInvTraitRateMatrix;

        if(matrixParams) {
            double[] tempInvVarianceMatrix = invVarianceMatrix;
            invVarianceMatrix = storedInvVarianceMatrix;
            storedInvVarianceMatrix = tempInvVarianceMatrix;

            double[] tempVarianceMatrix = varianceMatrix;
            varianceMatrix = storedVarianceMatrix;
            storedVarianceMatrix = tempVarianceMatrix;
        }

        double[] tempRootValuesVec = rootValuesVec;
        rootValuesVec= storedRootValuesVec;
        storedRootValuesVec = tempRootValuesVec;
    }



}

