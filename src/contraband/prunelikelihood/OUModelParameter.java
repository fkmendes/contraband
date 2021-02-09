package contraband.prunelikelihood;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import contraband.math.LUDecompositionForArray;
import contraband.math.MatrixUtilsContra;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;


/*
 * This class contains matrix and arrays for OU model, JOU model and DOU model
 */
public class OUModelParameter extends CalculationNode {
    final public Input<MorphologicalData> traitInput = new Input<>("trait","Morphological data set.");
    // OU model
    final public Input<RealParameter> alphaInput = new Input<>("alpha","An array of (nTraits * nTraits) elements, representing selection strength.", Input.Validate.REQUIRED);
    final public Input<RealParameter> thetaInput = new Input<>("theta","An array of nTraits elements, representing optimum trait values, elements in Theta vector.", Input.Validate.REQUIRED);
    final public Input<Integer> optNrInput = new Input<>("optNr","Number of theta (vectors).");
    final public Input<IntegerParameter> optAssignInput = new Input<>("optAssign", "the opt assignment for each node in the tree.");
    // JOU model
    final public Input<SigmaMatrix> jumpVCVMatInput = new Input<>("jumpVCVMat","Variance-covariance matrix of the normal jump distribution.");
    final public Input<RealParameter> jumpMeanValuesInput = new Input<>("jumpMean", "jump mean vector");
    final public Input<IntegerParameter> jumpIndicatorsInput = new Input<>("jumpIndicator","An array of boolean parameters corresponding to each branch in the tree. TRUE \n" +
            "indicates that a jump takes place at the beginning of the branch.");


    // DOU model
    final public Input<RealParameter> dAlphaInput = new Input<>("dAlpha","An array of (nTraits * nTraits) elements in decorrelation rate matrix.");


    private double[] alphaMat;
    private RealMatrix alphaRM;

    private double[] jumpVCVMat;
    private double[] jumpMeanVec;
    private double[] thetaVec;
    private double[] dAlphaMat;

    private int nTraits;
    private int nSpecies;
    private int nodeNr;
    private int optNr;
    private double[] thetaVecList;
    private Integer[] jumpIndicators;
    private Integer[] optAssignment;

    private double[] pMat;
    private double[] inverseP;

    // for doing LUDecompostion
    private double[] lu;
    private int[] pivot;
    private boolean[] evenSingular;
    private double [] identityMatrix;
    private boolean singularMatrix;

    // for store and restore



    @Override
    public void initAndValidate() {
        MorphologicalData data = traitInput.get();
        nSpecies = data.getSpeciesNr();
        nTraits = data.getTotalTraitNr();
        // number of nodes in the species tree including root
        nodeNr = 2 * nSpecies - 1;

        // alpha matrix is specified by nTraits * nTraits elements
        alphaMat = new double[nTraits * nTraits];
        alphaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        if (alphaInput.get().getDimension() != (nTraits * nTraits)){
            alphaInput.get().setDimension(nTraits * nTraits);
        }
        alphaMat = alphaInput.get().getDoubleValues();


        // there are optNr theta vectors in the tree
        // corresponding to categories assignment of each branch
        optNr = optNrInput.get();
        optAssignment = new Integer[nodeNr];
        if(optAssignInput.get().getDimension() != nodeNr) {
            optAssignInput.get().setDimension(nodeNr);
        }
        optAssignment = optAssignInput.get().getValues();

        // a list of theta vectors, including optNr different theta vectors
        // the index corresponds to the rate category
        thetaVecList = new double[nTraits * optNr];
        thetaVec = new double[nTraits];
        if(thetaInput.get().getDimension() != (nTraits * optNr)) {
            thetaInput.get().setDimension(nTraits * optNr);
        }
        // thetaVecList is a matrix with optNr row and nTraits column
        thetaVecList = thetaInput.get().getDoubleValues();

        // Initiate and validate inputs for JOU model
        if(jumpVCVMatInput.get() != null) {
            initJOUModel();
        }

        // Initiate and validate inputs for DOU model
        if(dAlphaInput.get() != null){
            initDOUModel();
        }

        // Initiate for Decomposition
        initLUDecomposition();
    }

    // getters
    public double[] getAlphaMat() { return  alphaMat; }

    public double[] getDAlphaMat() { return dAlphaMat; }

    public double[] getJumpVCVMat() { return jumpVCVMat; }

    public double[] getThetaVecForNode(int nodeIdx) {
        int index = optAssignment[nodeIdx];
        MatrixUtilsContra.getMatrixRow(thetaVecList, index, nTraits, thetaVec);
        return thetaVec;
    }

    public int getJumpIndicatorForNode(int nodeIdx) { return jumpIndicators[nodeIdx]; }

    public double[] getJumpMeanVec() {return jumpMeanVec;}

    public double[] getAlphaDecomposePMat() { return pMat; }

    public double[] getAlphaDecomposeInversePMat() { return inverseP; }

    /*
     * Initiate and validate inputs for JOU model
     */
    public void initJOUModel(){
        jumpVCVMat = jumpVCVMatInput.get().getSigmaMatrix();

        jumpMeanVec = new double[nTraits];
        if(jumpMeanValuesInput.get().getDimension() != nTraits) {
            jumpMeanValuesInput.get().setDimension(nTraits);
        }
        jumpMeanVec = jumpMeanValuesInput.get().getDoubleValues();

        jumpIndicators = new Integer[nodeNr];
        if(jumpIndicatorsInput.get().getDimension() != nodeNr) {
            jumpIndicatorsInput.get().setDimension(nodeNr);
        }
        jumpIndicators = jumpIndicatorsInput.get().getValues();
    }

    /*
     * Initiate and validate inputs for DOU model
     */
    public void initDOUModel(){
        dAlphaMat = new double[nTraits * nTraits];
        if (dAlphaInput.get().getDimension() != (nTraits * nTraits)){
            dAlphaInput.get().setDimension(nTraits * nTraits);
        }
        dAlphaMat = dAlphaInput.get().getDoubleValues();
    }

    /*
     * Initiate for LUDecomposition
     */
    private void initLUDecomposition() {
        lu = new double [nTraits * nTraits];
        pivot = new int[nTraits];
        evenSingular = new boolean[2];
        identityMatrix = new double[nTraits * nTraits];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, nTraits);
        }

        pMat = new double[nTraits * nTraits];
        inverseP = new double[nTraits * nTraits];
    }

    public void performAlphaEigenDecomposition() {
        if(dAlphaInput.get() != null){
            populateAlphaMatrix(alphaRM, dAlphaMat, nTraits);
        } else {
            populateAlphaMatrix(alphaRM, alphaMat, nTraits);
        }

        EigenDecomposition decompositionH = new EigenDecomposition(alphaRM);
        // normalize eigen vectors
        // NOTE: in R, eigen() returns vectors in decreasing order according to eigen values
        for (int i = 0; i < nTraits; i ++) {
            RealVector v = decompositionH.getEigenvector(i);
            double sum = 0.0;
            for (int j = 0; j < nTraits; j++){
                sum += (v.getEntry(j) * v.getEntry(j));
            }

            //pMat.setColumnVector(i,v.mapMultiply(1.0/ FastMath.sqrt(sum)));
            for(int k = 0; k < nTraits; k++){
                double value = v.getEntry(k) * (1.0/ FastMath.sqrt(sum));
                MatrixUtilsContra.setMatrixEntry(pMat, k, i, value, nTraits);
            }
        }

        // calculate the inverse matrix of pMat
        LUDecompositionForArray.ArrayLUDecomposition(pMat, lu, pivot, evenSingular, nTraits);
        try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evenSingular[1], nTraits, inverseP);
        }
        catch (SingularMatrixException e) {
            singularMatrix = true;
        }
    }

    public void updateOUModelParams() {
        if(alphaInput.isDirty()){
            alphaMat = alphaInput.get().getDoubleValues();
        }

        if(optAssignInput.isDirty()){
            optAssignment = optAssignInput.get().getValues();
        }

        if(thetaInput.isDirty()){
            thetaVecList = thetaInput.get().getDoubleValues();
        }

        if(jumpVCVMatInput.isDirty()){
            jumpVCVMatInput.get().updateParameters(false);
            jumpVCVMat = jumpVCVMatInput.get().getSigmaMatrix();
        }

        if(jumpIndicatorsInput.isDirty()){
            jumpIndicators = jumpIndicatorsInput.get().getValues();
        }

        if(jumpMeanValuesInput.isDirty()) {
            jumpMeanVec = jumpMeanValuesInput.get().getDoubleValues();
        }

        if(dAlphaInput.isDirty()){
            dAlphaMat = dAlphaInput.get().getDoubleValues();
        }
    }

    @Override
    public void store() {
        super.store();

    }

    @Override
    public void restore() {
        super.restore();

    }


    public static void populateAlphaMatrix(RealMatrix rm, double[] values, int nTraits) {
        for (int i = 0; i < nTraits; i++) {
            for (int j = 0; j < nTraits; j++) {
                rm.setEntry(i, j, MatrixUtilsContra.getMatrixEntry(values, i, j, nTraits));
            }
        }
    }
}
