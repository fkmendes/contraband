package contraband.math;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import contraband.prunelikelihood.MorphologicalData;
import contraband.prunelikelihood.SigmaMatrix;
import org.apache.commons.math3.util.FastMath;

public class GeneralNodeMath extends CalculationNode {
    final public Input<MorphologicalData> traitInput = new Input<>("trait","Morphological data set.");
    final public Input<SigmaMatrix> rateMatrixInput = new Input<>("rateMatrix", "Trait evolutionary variance-covariance matrix, i.e. trait rate matrix.", Input.Validate.REQUIRED);
    final public Input<SigmaMatrix> popMatrixInput = new Input<>("popMatrix", "Intraspecific variance-covariance matrix, i.e. population noise.");
    final public Input<Tree> treeInput = new Input<>("tree", "Phylogenetic tree.");
    final public Input<Boolean> shareRhoInput = new Input<>("shareCorrelation", "TRUE, if correlations are shared between rateMatrix and popMatrix.", true);
    final public Input<RealParameter> rootValuesInput = new Input<>("rootValues", "Trait values at the root.");


    private Tree tree;
    private int nodeNr;

    private boolean shareRho;
    private boolean matrixParams;
    private boolean sampleRoot;

    private SigmaMatrix rateMatrix;
    private SigmaMatrix popMatrix;

    private int nTraits;
    private int matDim;

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

    // this is for trees with sampled ancestors
    private double likForSA;
    private double[] traitVecForSA;

    private double [] traitRateMatrix;
    private double detTraitRateMat;
    private double [] invTraitRateMatrix;
    private double detInvTraitRateMat;

    private double storedDetTraitRateMat;
    private double [] storedInvTraitRateMatrix;
    private double storedDetInvTraitRateMat;
    private double [] storedTraitRateMatrix;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        nodeNr = tree.getNodeCount();

        nTraits = traitInput.get().getTotalTraitNr();
        matDim = nTraits * nTraits;

        // if trait rate matrix and intraspecific matrix share correlation
        shareRho = shareRhoInput.get();

        // get trait rate matrix
        rateMatrix = rateMatrixInput.get();

        // get intraspecific variance-covariance matrix
        if(popMatrixInput.get() == null){
            // ignore intraspecific variance
            matrixParams = false;
        } else {
            // account for intraspecific variance
            popMatrix = popMatrixInput.get();
            matrixParams = !rateMatrix.getOneRateOnly() || !popMatrix.getOneRateOnly();
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

        if(matrixParams) {
            // A C E are matrix for each node
            aMat = new double[matDim];
            cMat = new double[matDim];
            eMat = new double[matDim];
            aArray = new double[matDim * nodeNr];
            cArray = new double[matDim * nodeNr];
            eArray = new double[matDim * nodeNr];
        } else {
            // A C E are single double values for each node
            aArray = new double[nodeNr];
            cArray = new double[nodeNr];
            eArray = new double[nodeNr];
        }
    }

    private void initLmrArray() {
        // each nodes has a vector m
        mVecArray = new double[nTraits * nodeNr];
        mVec = new double[nTraits];
        mVecInit = new double[nTraits];
        // each nodes has a single double r
        rArray = new double[nodeNr];

        if(matrixParams) {
            // each node has a L matrix
            lMat = new double[matDim];
            lMatInit = new double[matDim];
            lArray = new double[matDim * nodeNr];
        } else {
            // each node has a single double value L
            lArray = new double[nodeNr];
        }
    }

    private void initUDecomposition() {
        lu = new double [matDim];
        pivot = new int[nTraits];
        evenSingular = new boolean[2];
        identityMatrix = new double[matDim];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, nTraits);
        }
    }

    private void initVarianceAndExpectation() {
        // each node has its variance = branch rate * branch length * trait rate matrix + intraspecific variance (if necessary)
        if(matrixParams){
            nodeVariance = new double[nodeNr * matDim];
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
        traitRateMatrix = new double[matDim];
        invTraitRateMatrix = new double[matDim];
        storedTraitRateMatrix = new double[matDim];
        storedInvTraitRateMatrix = new double[matDim];
    }

    //getters
    public double getAForNode (int nodeIdx) { return aArray[nodeIdx]; }

    public double[] getAMatForNode (int nodeIdx) {
        MatrixUtilsContra.getMatrixRow(aArray, nodeIdx, matDim, aMat);
        return aMat;
    }

    public double getCForNode (int nodeIdx) { return cArray[nodeIdx]; }

    public double[] getCMatForNode (int nodeIdx) {
        MatrixUtilsContra.getMatrixRow(cArray, nodeIdx, matDim, cMat);
        return cMat;
    }

    public double getEForNode (int nodeIdx) { return eArray[nodeIdx]; }

    public double[] getEMatForNode (int nodeIdx) {
        MatrixUtilsContra.getMatrixRow(eArray, nodeIdx, matDim, eMat);
        return eMat;
    }

    public double getfForNode (int nodeIdx) { return fArray[nodeIdx]; }

    public double getLForNode (int nodeIdx) { return lArray[nodeIdx]; }

    public double[] getLMatForNode (int nodeIdx) {
        MatrixUtilsContra.getMatrixRow(lArray, nodeIdx, matDim, lMat);
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

    public double[] getTraitsVec () {return traitsVec; }

    public double[] getTempVec () { return mVec; }

    public double[] getTraitRateMatrix () { return traitRateMatrix; }

    public double getTraitRateMatrixInverseDeterminant () { return detInvTraitRateMat; }

    // setters
    public void setLikelihoodForSampledAncestors(double value) {
        likForSA = value;
    }

    public void setAForNode (int nodeIdx, double value) { aArray[nodeIdx] = value; }

    public void setCForNode (int nodeIdx, double value) { cArray[nodeIdx] = value; }

    public void setEForNode (int nodeIdx, double value) { eArray[nodeIdx] = value; }

    public void setfForNode (int nodeIdx, double value) { fArray[nodeIdx] = value; }

    public void setLForNode (int nodeIdx, double value) { lArray[nodeIdx] = value; }

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
    public void operateOnTraitRateMatrix() {
        singularMatrix = false;

        traitRateMatrix = rateMatrix.getSigmaMatrix();

        // LUDecomposition
        LUDecompositionForArray.ArrayLUDecomposition(traitRateMatrix, lu, pivot, evenSingular, nTraits);

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
        System.arraycopy(invTraitRateMatrix, 0, storedInvTraitRateMatrix, 0, matDim);
        System.arraycopy(traitRateMatrix, 0, storedTraitRateMatrix, 0, matDim);
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

        double[] tempTraitRateMatrix = traitRateMatrix;
        traitRateMatrix = storedTraitRateMatrix;
        storedTraitRateMatrix = tempTraitRateMatrix;

        double[] tempInvTraitRateMatrix = invTraitRateMatrix;
        invTraitRateMatrix = storedInvTraitRateMatrix;
        storedInvTraitRateMatrix = tempInvTraitRateMatrix;

        double[] tempRootValuesVec = rootValuesVec;
        rootValuesVec= storedRootValuesVec;
        storedRootValuesVec = tempRootValuesVec;
    }



}

