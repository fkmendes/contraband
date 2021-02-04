package contraband.prunelikelihood;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.math.distributions.WishartDistribution;
import contraband.math.LUDecompositionForArray;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;

public class InverseMatrixGibbsOperator extends Operator {
    final public Input<Integer> dimensionInput = new Input<>("dimension", "The number of rows (columns) of the inverse matrix, representing the number of characters.", Input.Validate.REQUIRED);
    final public Input<RealParameter> valuesInput = new Input<>("values", "An array of values in the inverse matrix, representing either trait correlation or variance-covariance.", Input.Validate.REQUIRED);
    final public Input<WishartDistribution> wishartDistributionInput = new Input<>("distr", "The wishart distribution to randomly draw samples from." , Input.Validate.REQUIRED);
    final public Input<LiabilityLikelihood> likelihoodInput = new Input<>("likelihood", "Multivariate normal distribution for liabilities" , Input.Validate.REQUIRED);
    final public Input<NodeMath> nodeMathInput = new Input<>("nodeMath", "Multivariate normal distribution for liabilities" , Input.Validate.REQUIRED);

    private WishartDistribution wishartDistribution;
    private int numberOfTraits;
    private double[][] proposal;
    private int numberOfSpecies;
    private LiabilityLikelihood likelihood;
    private double [] priorInverseScaleMatrix;

    private double priorDf;
    private double[][] scaleMatrix;

    // for doing LUDecompostion
    private double[] lu;
    private int[] pivot;
    private boolean[] evensingular;
    private double [] identityMatrix;
    private boolean singularMatrix;

    private double[] S;
    private double[] S2;
    private double[] S2Plus;
    private double[] inverseS2Plus;
    private double[] data;

    @Override
    public void initAndValidate() {
        wishartDistribution = wishartDistributionInput.get();
        numberOfTraits = dimensionInput.get();
        if(wishartDistribution.df.get() < numberOfTraits) {
            throw new RuntimeException("InverseMatrixGibbsOperator:: Small n for large df.");
        }
        valuesInput.get().setDimension(numberOfTraits * numberOfTraits);
        proposal = new double[numberOfTraits][numberOfTraits];
        likelihood = likelihoodInput.get();
        priorDf = wishartDistribution.df.get();

        priorInverseScaleMatrix = new double[numberOfTraits * numberOfTraits];
        lu = new double [numberOfTraits * numberOfTraits];
        pivot = new int[numberOfTraits];
        evensingular = new boolean[2];
        identityMatrix = new double[numberOfTraits * numberOfTraits];

        // create an identity matrix for LUDecomposition
        for (int i = 0; i < numberOfTraits; i++) {
            MatrixUtilsContra.setMatrixEntry(identityMatrix, i, i, 1.0, numberOfTraits);
        }

        double[] priorScaleMatrix = wishartDistribution.scaleMatrix.get().getDoubleValues();
        LUDecompositionForArray.ArrayLUDecomposition(priorScaleMatrix, lu, pivot, evensingular, numberOfTraits);
        try {
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evensingular[1], numberOfTraits, priorInverseScaleMatrix);
        } catch (RuntimeException e) {
            singularMatrix = true;
        }

        // calculate sum-of-the-weighted-squares matrix over tree
        S = new double[numberOfTraits * numberOfTraits];
        numberOfSpecies = nodeMathInput.get().getNSpecies();
        S2 = new double[numberOfTraits * numberOfTraits];
        S2Plus = new double[numberOfTraits * numberOfTraits];
        inverseS2Plus = new double[numberOfTraits * numberOfTraits];
        data = new double[numberOfTraits];

        scaleMatrix = new double[numberOfTraits][numberOfTraits];
    }


    @Override
    public double proposal() {
        double[] scaleMatrixArr = getOperationScaleMatrixAndSetObservationCount();

        for (int i = 0; i < numberOfTraits; i++) {
            for(int j = 0; j < numberOfTraits; j++) {
                scaleMatrix[i][j] = MatrixUtilsContra.getMatrixEntry(scaleMatrixArr, i, j, numberOfTraits);
            }
        }

        proposal = WishartDistribution.nextWishart(priorDf, scaleMatrix);
        for(int i = 0; i < numberOfTraits; i ++){
            for (int j = 0; j < numberOfTraits; j ++){
                valuesInput.get().setValue(i * numberOfTraits + j, proposal[i][j]);
            }
        }

        if(singularMatrix){
            return Double.NEGATIVE_INFINITY;
        } else {
            return 0.0;
        }
    }

    private double[] getOperationScaleMatrixAndSetObservationCount() {

        // is a normal-normal-wishart model
        incrementOuterProduct(S, likelihood);

        try {

            System.arraycopy(S, 0, S2, 0, numberOfTraits * numberOfTraits);

            //S2 = priorInverseScaleMatrix.add(S2);
            MatrixUtilsContra.vectorAdd(S2, priorInverseScaleMatrix, S2Plus);

            //inverseS2 = (SymmetricMatrix) S2.inverse();
            LUDecompositionForArray.ArrayLUDecomposition(S2Plus, lu, pivot, evensingular, numberOfTraits);
            LUDecompositionForArray.populateInverseMatrix(lu, pivot, identityMatrix, evensingular[1], numberOfTraits, inverseS2Plus);
        } catch (RuntimeException e) {
            singularMatrix = true;
        }

        return inverseS2Plus;
    }

    private void incrementOuterProduct(double[]S, LiabilityLikelihood likelihood) {
        // the root values
        // likelihood.getDistribution().getMean();
        double[] mean = nodeMathInput.get().getRootValuesArr();


        // the traitValues;
        //List<Attribute<double[]>> dataList = likelihood.getDataList();
        double[] traitValues = likelihood.getCombinedTraitDataArr().clone();

        // iterate each species
        for (int i = 0; i < numberOfSpecies; i++) {
            // trait values for i-th species
             MatrixUtilsContra.getMatrixRow(traitValues, i, numberOfTraits, data);

             // subtract mean values
            for (int j = 0; j < numberOfTraits; j++) {
                data[j] -= mean[j];
            }

            // accumulate
            for (int k = 0; k < numberOfTraits; k++) {
                for (int l = k; l < numberOfTraits; l++) {
                    //S[l][k] = S[k][l] += data[k] * data[l];
                    S[k * numberOfTraits + l] = S[l * numberOfTraits + k] +=  data[k] * data[l];
                }
            }
        }
    }
}
