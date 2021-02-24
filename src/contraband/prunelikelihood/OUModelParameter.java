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
    private SigmaMatrix jumpMatrix;
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
    private double[] pMatTranspose;
    private double[] inversePTranspose;
    private double[] eigenValues;

    // for doing LUDecompostion
    private double[] lu;
    private int[] pivot;
    private boolean[] evenSingular;
    private double [] identityMatrix;
    private boolean singularMatrix;

    // for store and restore
    private double[] storedThetaVecList;
    private Integer[] storedOptAssignment;
    private double[] storedAlphaMat;
    private double[] storedPMat;
    private double[] storedInverseP;
    private double[] storedPMatTranspose;
    private double[] storedInversePTranspose;
    private double[] storedEigenValues;
    private boolean storedSingular;
    private double[] storedJumpVCVMat;
    private Integer[] storedJumpIndicators;
    private double[] storedJumpMeanVec;
    private double[] storedDAlphaMat;

    private String modelType;


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

        modelType = "OU";

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
            modelType = "JOU";
            initJOUModel();
        }

        // Initiate and validate inputs for DOU model
        if(dAlphaInput.get() != null){
            modelType = "DOU";
            initDOUModel();
        }

        // Initiate arrays and matrices
        initLUDecomposition();
        populateStoreAndRestore();
    }

    // getters
    public double[] getAlphaMat() { return  alphaMat; }

    public double[] getDAlphaMat() { return dAlphaMat; }

    public double[] getJumpVCVMat() { return jumpMatrix.getSigmaMatrix(); }

    public double[] getThetaVecForNode(int nodeIdx) {
        int index = optAssignment[nodeIdx];
        MatrixUtilsContra.getMatrixRow(thetaVecList, index, nTraits, thetaVec);
        return thetaVec;
    }

    public int getJumpIndicatorForNode(int nodeIdx) { return jumpIndicators[nodeIdx]; }

    public double[] getJumpMeanVec() {return jumpMeanVec;}

    public double[] getAlphaDecomposePMat() { return pMat; }

    public double[] getAlphaDecomposeInversePMat() { return inverseP; }

    public double[] getPMatTranspose() { return pMatTranspose; }

    public double[] getInversePTranspose() {return inversePTranspose; }

    public double[] getEigenValues() { return eigenValues; }

    public boolean isSingularAlphaMatrix() { return singularMatrix; }

    public String getModelType() { return modelType; }

    /*
     * Initiate and validate inputs for JOU model
     */
    public void initJOUModel(){
        jumpMatrix = jumpVCVMatInput.get();
        jumpVCVMatInput.get().populateSigmaMatrix();
        //jumpVCVMat = jumpVCVMatInput.get().getSigmaMatrix();

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
        eigenValues = new double[nTraits];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < nTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, nTraits);
        }

        pMat = new double[nTraits * nTraits];
        inverseP = new double[nTraits * nTraits];
        pMatTranspose = new double[nTraits * nTraits];
        inversePTranspose = new double[nTraits * nTraits];
    }

    // initiate for store and restore
    private void populateStoreAndRestore() {
        storedThetaVecList = new double[nTraits * optNr];
        storedOptAssignment = new Integer[nodeNr];
        storedAlphaMat = new double[nTraits * nTraits];
        storedPMat = new double[nTraits * nTraits];
        storedInverseP = new double[nTraits * nTraits];
        storedPMatTranspose = new double[nTraits * nTraits];
        storedInversePTranspose = new double[nTraits * nTraits];
        storedEigenValues = new double[nTraits];
        if(modelType.equals("JOU")){
            //storedJumpVCVMat = new double [jumpVCVMat.length];
            storedJumpIndicators = new Integer[nodeNr];
            storedJumpMeanVec = new double[nTraits];
        }
        if(modelType.equals("DOU")){
            storedDAlphaMat = new double[nTraits * nTraits];
        }
    }

    public void performAlphaEigenDecomposition() {
        singularMatrix = false;
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

        // get eigen values
        eigenValues = decompositionH.getRealEigenvalues();

        // calculate the inverse matrix of pMat
        LUDecompositionForArray.ArrayLUDecomposition(pMat, lu, pivot, evenSingular, nTraits);
        try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evenSingular[1], nTraits, inverseP);
        }
        catch (RuntimeException e) {
            singularMatrix = true;
        }

        // get the transpose matrix
        populateTransposeMatrix();
    }

    private void populateTransposeMatrix(){
        MatrixUtilsContra.matrixTranspose(pMat, nTraits, pMatTranspose);
        MatrixUtilsContra.matrixTranspose(inverseP, nTraits, inversePTranspose);
    }

    public void updateOUModelParams() {
        if(alphaInput.isDirty()){
            alphaMat = alphaInput.get().getDoubleValues();
            performAlphaEigenDecomposition();
        }

        if(optAssignInput.isDirty()){
            optAssignment = optAssignInput.get().getValues();
        }

        if(thetaInput.isDirty()){
            thetaVecList = thetaInput.get().getDoubleValues();
        }

        if(jumpIndicatorsInput.isDirty()){
            jumpIndicators = jumpIndicatorsInput.get().getValues();
        }

            if(jumpVCVMatInput.isDirty()){
            jumpMatrix.updateParameters(false);
            //jumpVCVMat = jumpMatrix.getSigmaMatrix();
           }

            if(jumpMeanValuesInput.isDirty()) {
                jumpMeanVec = jumpMeanValuesInput.get().getDoubleValues();
           }


        if(dAlphaInput.isDirty()){
            dAlphaMat = dAlphaInput.get().getDoubleValues();
            performAlphaEigenDecomposition();
        }
    }

    @Override
    public void store() {
        super.store();
        System.arraycopy(alphaMat, 0, storedAlphaMat, 0, alphaMat.length);

        System.arraycopy(eigenValues, 0, storedEigenValues , 0, eigenValues.length);

        System.arraycopy(pMat, 0, storedPMat, 0, pMat.length);

        System.arraycopy(pMatTranspose, 0, storedPMatTranspose, 0, pMatTranspose.length);

        System.arraycopy(inverseP, 0, storedInverseP, 0, inverseP.length);

        System.arraycopy(inversePTranspose, 0, storedInversePTranspose, 0, inversePTranspose.length);

        System.arraycopy(optAssignment, 0, storedOptAssignment, 0, optAssignment.length);

        System.arraycopy(thetaVecList, 0, storedThetaVecList, 0, thetaVecList.length);

        if(modelType.equals("JOU")){
            //System.arraycopy(jumpVCVMat, 0, storedJumpVCVMat, 0, jumpVCVMat.length);
            System.arraycopy(jumpIndicators, 0, storedJumpIndicators, 0, jumpIndicators.length);
            System.arraycopy(jumpMeanVec, 0, storedJumpMeanVec, 0, jumpMeanVec.length);
        }
        if(modelType.equals("DOU")){
            System.arraycopy(dAlphaMat, 0, storedDAlphaMat, 0, dAlphaMat.length);
        }
        storedSingular = singularMatrix;
    }

    @Override
    public void restore() {
        super.restore();

        double[] tempAlphaMat = alphaMat;
        alphaMat = storedAlphaMat;
        storedAlphaMat = tempAlphaMat;

        double[] tempEigenValues = eigenValues;
        eigenValues = storedEigenValues;
        storedEigenValues = tempEigenValues;

        double[] tempPMat = pMat;
        pMat = storedPMat;
        storedPMat = tempPMat;

        double[] tempPMatTranspose = pMatTranspose;
        pMatTranspose = storedPMatTranspose;
        storedPMatTranspose = tempPMatTranspose;

        double[] tempInverseP = inverseP;
        inverseP = storedInverseP;
        storedInverseP = tempInverseP;

        double[] tempInversePTranspose = inversePTranspose;
        inversePTranspose = storedInversePTranspose;
        storedInversePTranspose = tempInversePTranspose;

        Integer[] tempOptAssignment = optAssignment;
        optAssignment = storedOptAssignment;
        storedOptAssignment = tempOptAssignment;

        double[] tempThetaVecList = thetaVecList;
        thetaVecList = storedThetaVecList;
        storedThetaVecList = tempThetaVecList;

        boolean tempSingular = singularMatrix;
        singularMatrix = storedSingular;
        storedSingular = tempSingular;

        if(modelType.equals("JOU")){
            Integer[] tempJumpIndicators = jumpIndicators;
            jumpIndicators = storedJumpIndicators;
            storedJumpIndicators = tempJumpIndicators;

            //double[] tempJumpVCVMat = jumpVCVMat;
            //jumpVCVMat = storedJumpVCVMat;
            //storedJumpVCVMat = tempJumpVCVMat;

            double[] tempJumpMeanVec = jumpMeanVec;
            jumpMeanVec = storedJumpMeanVec;
            storedJumpMeanVec = tempJumpMeanVec;
        }
        if(modelType.equals("DOU")){
            double[] tempDAlphaMat = dAlphaMat;
            dAlphaMat = storedDAlphaMat;
            storedDAlphaMat = tempDAlphaMat;
        }
    }


    public static void populateAlphaMatrix(RealMatrix rm, double[] values, int nTraits) {
        for (int i = 0; i < nTraits; i++) {
            for (int j = 0; j < nTraits; j++) {
                rm.setEntry(i, j, MatrixUtilsContra.getMatrixEntry(values, i, j, nTraits));
            }
        }
    }
}
