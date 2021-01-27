package contraband.prunelikelihood;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.evolution.tree.Tree;
import beast.math.matrixalgebra.CholeskyDecomposition;
import beast.math.matrixalgebra.IllegalDimension;
import contraband.math.LUDecompositionForArray;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import outercore.parameter.KeyRealParameter;

// This class deals with different kinds of morphological traits, including continuous, binary discrete, ordered discrete and unordered discrete,
// and populates a combined data set using continuous trait values and liabilities of discrete traits.
public class MorphologicalData extends CalculationNode{
    final public Input<KeyRealParameter> traitsValuesInput = new Input<>("traits","Trait values at tips.");
    final public Input<BinaryDiscreteTraits> binaryDiscreteTraitsInput = new Input<>("binaryDiscreteTraits", "Object for binary discrete traits.");
    final public Input<OrderedDiscreteTraits> orderedDiscreteTraitsInput = new Input<>("orderedDiscreteTraits", "Object for ordered discrete traits.");
    final public Input<UnorderedDiscreteTraits> unorderedDiscreteTraitsInput = new Input<>("unorderedDiscreteTraits", "Object for unordered discrete traits.");
    final public Input<Tree> treeInput = new Input<>("tree", "Phylogenetic tree.");
    final public Input<Boolean> transformInput = new Input<>("transform", "TRUE, if data needs to be transformed to be independent", false);
    final public Input<NodeMath> nodeMathInput = new Input<>("nodeMath","Node information that will be used in PCM likelihood calculation.", Input.Validate.REQUIRED);

    private Tree  tree;
    private int totalTraitNr;
    private int speciesNr;

    private double[] traitValuesArr;
    private double[] storedTraitValuesArr;

    private int binaryTraitIndex;
    private int orderedTraitIndex;
    private  int unorderedLiabilityIndex;

    private int contTraitNr = 0;
    private int binaryTraitNr = 0;
    private int orderedTraitNr = 0;
    private int unorderedLiabilityNr = 0;

    private int matrixDimension;
    private boolean transformedData;
    private double[][] cholesky;
    private double[][] lowerMatrix;
    private double[] lowerMatArr;
    private double[] dataTransformMat;
    private double[] upperMatrix;
    private double[] transformedTraitValues;
    private double[] storedTransformedTraitValues;

    private double[] lu;
    private int[] pivot;
    private boolean[] evensingular;
    private double [] identityMatrix;
    private boolean singularMatrix;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        int leafNodeCount = tree.getLeafNodeCount();

        // collect trait information
        totalTraitNr = 0;
        speciesNr = 0;
       if (traitsValuesInput.get() != null) {
           KeyRealParameter traitsValues = traitsValuesInput.get();
           totalTraitNr += traitsValues.getMinorDimension1();
           speciesNr = traitsValues.getMinorDimension2();
       }
       if(binaryDiscreteTraitsInput.get() != null) {
           totalTraitNr += binaryDiscreteTraitsInput.get().getLiabilityNr();
           speciesNr = binaryDiscreteTraitsInput.get().getSpeciesNr();
       }
       if(orderedDiscreteTraitsInput.get() != null) {
           totalTraitNr += orderedDiscreteTraitsInput.get().getLiabilityNr();
           speciesNr = orderedDiscreteTraitsInput.get().getSpeciesNr();
       }
       if(unorderedDiscreteTraitsInput.get() != null){
           totalTraitNr += unorderedDiscreteTraitsInput.get().getNrOfLiabilities();
           speciesNr = unorderedDiscreteTraitsInput.get().getSpeciesNr();
       }

       if(speciesNr != leafNodeCount) {
           throw new RuntimeException("MorphologicalData::Species number in the data set does not match the leaf node number in tree.");
       }

        transformedData = transformInput.get();

        initArrays();
        initLUDecomposition();

    }

    // getters
    public int getTotalTraitNr() { return totalTraitNr; }

    public int getSpeciesNr () { return speciesNr; }

    public double[] getMorphologicalData () {
        if(transformedData) {
            return traitValuesArr;
        } else {
            return transformedTraitValues;
        }
    }

    public boolean getSingularTransformMatrix () { return singularMatrix; }

    public boolean getTransformDataFlag () { return transformedData; }


    // initiate
    private void initArrays() {
        matrixDimension = totalTraitNr * totalTraitNr;
        cholesky = new double[totalTraitNr][totalTraitNr];
        lowerMatrix = new double[totalTraitNr][totalTraitNr];
        lowerMatArr = new double[matrixDimension ];
        dataTransformMat = new  double[matrixDimension];
        upperMatrix = new double[matrixDimension];
        traitValuesArr = new double[speciesNr * totalTraitNr];
        transformedTraitValues = new double[speciesNr * totalTraitNr];
    }

    private void initLUDecomposition() {
        lu = new double [matrixDimension];
        pivot = new int[totalTraitNr];
        evensingular = new boolean[2];
        identityMatrix = new double[matrixDimension];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < totalTraitNr; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, totalTraitNr);
        }
    }


    // initialize the combined morphological data set
    private void populateTraitData(){
        /*
         * populate combined trait values in an array
         */
        traitValuesArr = new double[speciesNr * totalTraitNr];
        storedTraitValuesArr = new double[speciesNr * totalTraitNr];

        int index = 0;
        // (1) populate continuous trait values
        if(traitsValuesInput.get() != null) {
            populateContTraitArr(traitsValuesInput.get(), tree, contTraitNr, totalTraitNr, traitValuesArr);
            index += contTraitNr;
        }
        binaryTraitIndex = index;

        // (2) populate liabilities for binary discrete trait values
        if(binaryDiscreteTraitsInput.get() != null){
            double[] binaryLiabilities = binaryDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(binaryLiabilities, traitValuesArr, index, binaryTraitNr, totalTraitNr, speciesNr);
            index += binaryTraitNr;
        }
        orderedTraitIndex = index;

        // (3) populate liabilities for ordered discrete trait values
        if(orderedDiscreteTraitsInput.get() != null) {
            double[] orderedLiabilities = orderedDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(orderedLiabilities , traitValuesArr, index, orderedTraitNr, totalTraitNr, speciesNr);
            index += orderedTraitNr;
        }
        unorderedLiabilityIndex = index;

        // (4) populate liabilities for unordered discrete trait values
        if(unorderedDiscreteTraitsInput.get() != null) {
            double[] unorderedLiabilities = unorderedDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(unorderedLiabilities, traitValuesArr, index, unorderedLiabilityNr, totalTraitNr, speciesNr);
        }

        // (5) transform the data if specified
        if(transformedData) {
            transformTraitValues(traitValuesArr);
        }
    }

    private void updateTraitValuesArr(boolean updateRateMatrix){
        boolean update = false;
        // (1) update liabilities for binary discrete trait values
        if(binaryDiscreteTraitsInput.get()!= null && binaryDiscreteTraitsInput.get().liabilitiesInput.isDirty()){
            double[] binaryLiabilities = binaryDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(binaryLiabilities, traitValuesArr, binaryTraitIndex, binaryTraitNr, totalTraitNr, speciesNr);
            update = true;
        }

        // (2) update liabilities for ordered discrete trait values
        if(orderedDiscreteTraitsInput.get() != null && orderedDiscreteTraitsInput.get().liabilitiesInput.isDirty()) {
            double[] orderedLiabilities = orderedDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(orderedLiabilities, traitValuesArr, orderedTraitIndex, orderedTraitNr, totalTraitNr, speciesNr);
            update = true;
        }

        // (3) update liabilities for unordered discrete trait values
        if(unorderedDiscreteTraitsInput.get() != null && unorderedDiscreteTraitsInput.get().liabilitiesInput.isDirty()){
            double[] unorderedLiabilities = unorderedDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(unorderedLiabilities, traitValuesArr, unorderedLiabilityIndex, unorderedLiabilityNr, totalTraitNr, speciesNr);
            update = true;
        }

        // if liability or trait rate matrix is changed, we have to update the transformed trait values
        if(transformedData) {
            if (update || updateRateMatrix) {
                // (4) transform the original trait values based on the trait rate matrix
                // so that the transformed trait value are independent
                transformTraitValues(traitValuesArr);
            }
        }
    }



    private void transformTraitValues (double[] liabilities) {
        // (1) copy trait rate matrix to a 2D array
        for(int i = 0; i < totalTraitNr; i++){
            System.arraycopy(nodeMathInput.get().getTraitRateMatrix(), i * totalTraitNr, cholesky[i], 0, totalTraitNr);
        }

        // (2) perform cholesky decomposition, R = L * L.transpose
        try {
            lowerMatrix = (new CholeskyDecomposition(cholesky)).getL();
            // caution: this returns the lower triangular form
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Numerical exception in WishartDistribution");
        }
        // copy the 2D lower matrix to a double array
        for(int i = 0; i < totalTraitNr; i++){
            System.arraycopy(lowerMatrix[i], 0, lowerMatArr, i * totalTraitNr, totalTraitNr);
        }

        // (3) get the inverse matrix, A = L.inverse
        LUDecompositionForArray.ArrayLUDecomposition(lowerMatArr, lu, pivot, evensingular, totalTraitNr);
        try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evensingular[1], totalTraitNr, dataTransformMat);
        } catch (RuntimeException e) {
            singularMatrix = true;
        }

        // (4) Z = M * A.transpose
        // the transpose the matrix
        MatrixUtilsContra.matrixTranspose(dataTransformMat, totalTraitNr, upperMatrix);
        MatrixUtilsContra.matrixMultiply(liabilities, upperMatrix, speciesNr, totalTraitNr, transformedTraitValues);
    }


    /*
     * The following method can be static in utils.
     */
    // add continuous trait values to the traitValuesArr
    private void populateContTraitArr(KeyRealParameter traitValues, Tree tree, int contTraitNr, int totalTraitNr, double[] traitValuesArr) {
        // according to node number of tips
        for (int i = 0; i < tree.getLeafNodeCount(); i ++) {
            // get all traits values for this species
            Double[] traitForSpecies = traitValues.getRowValues(tree.getNode(i).getID());
            for (int j= 0; j < contTraitNr; j ++) {
                // populate the traits one by one in an array
                traitValuesArr[i*totalTraitNr + j] = traitForSpecies[j];
            }
        }
    }

    // add liabilities of discrete traits to the traitValuesArr
    private void populateCombinedTraitValuesArr(double[] liabilities, double[] traitValuesArr, int index, int traitNr, int totalNr, int speciesNr){
        for (int i = 0; i < speciesNr; i++) {
            System.arraycopy(liabilities, i * traitNr, traitValuesArr, i * totalNr + index, traitNr);
        }
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    public void store() {
        super.store();
        if (transformedData){
            System.arraycopy(transformedTraitValues, 0, storedTransformedTraitValues, 0, speciesNr * totalTraitNr);
        } else {
            System.arraycopy(traitValuesArr, 0, storedTraitValuesArr, 0, speciesNr * totalTraitNr);
        }
    }

    @Override
    public void restore() {
        super.restore();
        if(transformedData) {
            double[] tempTraitValuesArray = transformedTraitValues;
            transformedTraitValues = storedTransformedTraitValues;
            storedTransformedTraitValues = tempTraitValuesArray;
        } else {
            double[] tempTraitValuesArr = traitValuesArr;
            traitValuesArr = storedTraitValuesArr;
            storedTraitValuesArr = tempTraitValuesArr;
        }
    }


}
