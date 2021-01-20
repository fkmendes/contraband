package contraband.prunelikelihood;

import beast.core.Input;
import beast.evolution.tree.Tree;
import beast.math.matrixalgebra.CholeskyDecomposition;
import beast.math.matrixalgebra.IllegalDimension;
import contraband.math.LUDecompositionForArray;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.utils.PruneLikelihoodUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;

public class TransformedLiabilityLikelihood extends PruneLikelihoodProcess {
    final public Input<BinaryDiscreteTraits> binaryDiscreteTraitsInput = new Input<>("binaryDiscreteTraits", "Object for binary discrete traits.");
    final public Input<OrderedDiscreteTraits> orderedDiscreteTraitsInput = new Input<>("orderedDiscreteTraits", "Object for ordered discrete traits.");
    final public Input<UnorderedDiscreteTraits> unorderedDiscreteTraitsInput = new Input<>("unorderedDiscreteTraits", "Object for unordered discrete traits.");

    private Tree tree;
    private int contTraitNr = 0;
    private int binaryTraitNr = 0;
    private int orderedTraitNr = 0;
    private int unorderedLiabilityNr = 0;
    private int totalTraitNr;

    private int nSpecies;
    private double[] traitValuesArr;


    private int binaryTraitIndex;
    private int orderedTraitIndex;
    private  int unorderedLiabilityIndex;

    RealMatrix rateRealMatrix;

    @Override
    public void initAndValidate() {
        // get the tree
        tree = treeInput.get();
        nSpecies = treeInput.get().getLeafNodeCount();

        // get number of liabilities for discrete traits
        if(binaryDiscreteTraitsInput.get() != null){
            binaryTraitNr = binaryDiscreteTraitsInput.get().getLiabilityNr();
        }
        if(orderedDiscreteTraitsInput.get() != null){
            orderedTraitNr = orderedDiscreteTraitsInput.get().getLiabilityNr();
        }
        if(unorderedDiscreteTraitsInput.get() != null){
            UnorderedDiscreteTraits unorderedDiscreteTraits = unorderedDiscreteTraitsInput.get();
            unorderedLiabilityNr = unorderedDiscreteTraits.getNrOfLiabilities();
        }
        // get number of continuous traits
        if(traitsValuesInput.get() != null) {
            contTraitNr = traitsValuesInput.get().getMinorDimension1();
        }

        // get total number of traits
        totalTraitNr = binaryTraitNr + orderedTraitNr + unorderedLiabilityNr + contTraitNr;
        if (totalTraitNr == 0) {
            throw new RuntimeException("LiabilityLikelihood::At least one kind of data should be input.");
        }
        setNTraits(totalTraitNr);
        rateRealMatrix = new Array2DRowRealMatrix(new double[totalTraitNr][totalTraitNr]);

        super.initAndValidate();

        getNodeMath().populateTraitRateMatrix();
        getNodeMath().performMatrixOperations();
        getNodeMath().populateTransformedTraitValues(traitValuesArr);
        setTraitValuesArr(getNodeMath().getTransformedTraitValues());
    }

    // add liabilities of discrete traits to the traitValuesArr
    private void populateCombinedTraitValuesArr(double[] liabilities, double[] traitValuesArr, int index, int traitNr, int totalNr, int speciesNr){
        for (int i = 0; i < speciesNr; i++) {
            System.arraycopy(liabilities, i * traitNr, traitValuesArr, i * totalNr + index, traitNr);
        }
    }

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

    @Override
    protected void populateTraitData(){
        /*
         * populate combined trait values in an array
         */
        traitValuesArr = new double[nSpecies * totalTraitNr];

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
            populateCombinedTraitValuesArr(binaryLiabilities, traitValuesArr, index, binaryTraitNr, totalTraitNr, nSpecies);
            index += binaryTraitNr;
        }
        orderedTraitIndex = index;

        // (3) populate liabilities for ordered discrete trait values
        if(orderedDiscreteTraitsInput.get() != null) {
            double[] orderedLiabilities = orderedDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(orderedLiabilities , traitValuesArr, index, orderedTraitNr, totalTraitNr, nSpecies);
            index += orderedTraitNr;
        }
        unorderedLiabilityIndex = index;

        // (4) populate liabilities for unordered discrete trait values
        if(unorderedDiscreteTraitsInput.get() != null){
            double[] unorderedLiabilities = unorderedDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(unorderedLiabilities, traitValuesArr, index, unorderedLiabilityNr, totalTraitNr, nSpecies);
        }

    }

    public double[] getCombinedTraitDataArr() { return traitValuesArr; }

    @Override
    public double calculateLogP() {
        boolean updateTraitRateMatrix = false;
        if(nodeMathInput.isDirty()) {
            //update parameters if some parameters are dirty
            // inverse rho matrix
            // sigma value(s)
            // any variables related to trait rate matrix
            updateTraitRateMatrix = getNodeMath().updateParameters();
        }
        if(updateTraitRateMatrix) {
            // populate trait rate matrix with updated parameters
            getNodeMath().populateTraitRateMatrix();

            // if the trait rate matrix is nearly singular
            // it should be rejected in advance
            if(checkNearlySingularForRateMatrix()) {return Double.NEGATIVE_INFINITY;}

            // four quantities will be required in likelihood calculation
            // trait rate matrix, inverse trait rate matrix
            // det(trait rate matrix), det(inverse trait rate matrix)
            getNodeMath().performMatrixOperations();
        }

        // update trait values if some liability is dirty
        updateTraitValuesArr(updateTraitRateMatrix);

        super.populateLogP();

        return getLogP();
    }

    private void updateTraitValuesArr(boolean updateRateMatrix){
        boolean update = false;
        // (1) update liabilities for binary discrete trait values
        if(binaryDiscreteTraitsInput.isDirty()){
            double[] binaryLiabilities = binaryDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(binaryLiabilities, traitValuesArr, binaryTraitIndex, binaryTraitNr, totalTraitNr, nSpecies);
            update = true;
        }

        // (2) update liabilities for ordered discrete trait values
        if(orderedDiscreteTraitsInput.isDirty()) {
            double[] orderedLiabilities = orderedDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(orderedLiabilities, traitValuesArr, orderedTraitIndex, orderedTraitNr, totalTraitNr, nSpecies);
            update = true;
        }

        // (3) update liabilities for unordered discrete trait values
        if(unorderedDiscreteTraitsInput.isDirty()){
            double[] unorderedLiabilities = unorderedDiscreteTraitsInput.get().getLiabilities();
            populateCombinedTraitValuesArr(unorderedLiabilities, traitValuesArr, unorderedLiabilityIndex, unorderedLiabilityNr, totalTraitNr, nSpecies);
            update = true;
        }

        // if liability or trait rate matrix is changed, we have to update the transformed trait values
        if(update || updateRateMatrix) {
            // (4) transform the original trait values based on the trait rate matrix
            // so that the transformed trait value are independent
            //populateTransformedTraitValues(traitValuesArr);
            getNodeMath().populateTransformedTraitValues(traitValuesArr);

            // (5) update the trait values in the PruneLikelihoodProcess class
            setTraitValuesArr(getNodeMath().getTransformedTraitValues());
        }
    }

    @Override
    protected void calculateLmrForTips(NodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        PruneLikelihoodUtils.populateLmrForTipWithShrinkage(nodeMath, traitValuesArr, nTraits, nodeIdx);
    }

    @Override
    protected void calculateLmrForInternalNodes(NodeMath nodeMath, int nTraits, int nodeIdx) {
        PruneLikelihoodUtils.populateLmrForInternalNodeWithShrinkage(nodeMath, nTraits, nodeIdx);
    }

    @Override
    protected double calculateLikelihood(NodeMath nodeMath, double l0, double[] m0, double r0, int rootIdx){
        //double vCD = nodeMath.getVarianceForNode(rootIdx);
        //double root2Subtract = -0.5 * nodeMath.getTraitRateMatrixDeterminant() - ((getNTraits() / 2.0) * Math.log(2 * Math.PI * vCD));
        return MatrixUtilsContra.vecTransScalarMultiply(nodeMath.getRootValuesArr(),
                l0, getNTraits()) +
                MatrixUtilsContra.vectorDotMultiply(nodeMath.getRootValuesArr(), m0) +
                r0 +
                nodeMath.getLikelihoodForSampledAncestors();

    }

    private boolean checkNearlySingularForRateMatrix () {
        // initialize
        boolean nearlySingularRateMatrix = false;

        // for one trait analysis
        if (getNTraits() == 1) {
            if (getNodeMath().getTraitRateMatrix()[0] < 1.0E-5) {
                nearlySingularRateMatrix = true;
            }
        }

        // for multiple traits
        else {
            for (int i = 0; i < getNTraits(); i++) {
                MatrixUtilsContra.getMatrixRow(getNodeMath().getTraitRateMatrix(), i, getNTraits(), getNodeMath().getRateMatrixRow());
                rateRealMatrix.setRow(i, getNodeMath().getRateMatrixRow());
            }


            double[] singularValues = new SingularValueDecomposition(rateRealMatrix).getSingularValues();
            double min = Arrays.stream(singularValues).min().getAsDouble();
            double max = Arrays.stream(singularValues).max().getAsDouble();

            double[] eValues = new EigenDecomposition(rateRealMatrix).getRealEigenvalues();

            for (double ei : eValues) {
                if (ei < 1.0E-5) {
                    nearlySingularRateMatrix = true;
                }
            }

            if ((min / max) < 1.0E-6) {
                nearlySingularRateMatrix = true;
            }
        }

        return nearlySingularRateMatrix;
    }
}
