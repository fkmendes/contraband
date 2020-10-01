package contraband.prunelikelihood;

import beast.core.Description;
import beast.core.Input;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.utils.PruneLikelihoodUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import java.util.Arrays;

@Description("This class implements likelihood for continuous traits under Brownian model.\n" +
        "The calculation uses Venelin's PCM likelihood.")


public class BMPruneLikelihood extends PruneLikelihoodProcess {
    final public Input<Boolean> includeRootInput = new Input<>("includeRoot", "TRUE, if the likelihood at will be subtracted.", false);

    private RealMatrix rateRealMatrix;
    private boolean nearlySingularRateMatrix;
    private boolean includeRoot;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // initialize a RealMatrix to check if the trait rate matrix is nearly singular or not
        rateRealMatrix = new Array2DRowRealMatrix(new double[getNTraits()][getNTraits()]);

        getNodeMath().populateTraitRateMatrix();
        getNodeMath().performMatrixOperations();

        includeRoot = includeRootInput.get();

        setPopSE(false);
    }

    @Override
    public double calculateLogP() {
        boolean updateTraitRateMatrix = false;
        if(nodeMathInput.isDirty()) {
            // update parameters if some parameters are dirty
            updateTraitRateMatrix = getNodeMath().updateParameters();
        }
        if(updateTraitRateMatrix) {
            // get the variance and covariance of traits
            // and populate a 1D double array of RateMatrix
            getNodeMath().populateTraitRateMatrix();

            // check if the populated RateMatrix is nearly singular or not
            // by creating a RealMatrix object
            checkNearlySingularForRateMatrix();

            // if the RateMatrix is nearly singular
            // reject this state
            if (nearlySingularRateMatrix) {
                return Double.NEGATIVE_INFINITY;
            }

            // populate the inverse matrix of the RateMatrix
            // and calculate the determinant of the RateMatrix
            getNodeMath().performMatrixOperations();

            // if the RateMatrix is singular
            // reject this state
            if (getNodeMath().isSingularMatrix()) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        super.populateLogP();

        // if TRUE, the likelihood at the root will be subtracted
        // so that the likelihood matches mcmcTree without using shrinkage method
        if(includeRoot) {
            double vCD = getNodeMath().getVarianceForNode(getRootIndex());
            double root2Subtract = -0.5 * getNodeMath().getTraitRateMatrixDeterminant() - ((getNTraits() / 2.0) * Math.log(2 * Math.PI * vCD));
            return getLogP() - root2Subtract;
        }

        return getLogP();
    }

    private void checkNearlySingularForRateMatrix () {
        // initialize
        nearlySingularRateMatrix = false;

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
    }


    @Override
    protected void calculateLmrForTips(NodeMath nodeMath, double[] traitValuesArr, int nTraits, int nodeIdx) {
        PruneLikelihoodUtils.populateLmrForTip(nodeMath, traitValuesArr, nTraits, nodeIdx);
    }

    @Override
    protected void calculateLmrForInternalNodes(NodeMath nodeMath, int nTraits, int nodeIdx) {
        PruneLikelihoodUtils.populateLmrForInternalNode(nodeMath, nTraits, nodeIdx);
    }

    @Override
    protected double calculateLikelihood(NodeMath nodeMath, double l0, double[] m0, double r0, int rootIdx){
        return  l0 * MatrixUtilsContra.tVecDotMatrixDotVec(
                nodeMath.getRootValuesArr(),
                nodeMath.getTraitRateMatrixInverse(),
                getNTraits()) +
                MatrixUtilsContra.vectorDotMultiply(
                        nodeMath.getRootValuesArr(),
                        m0) +
                r0 + nodeMath.getLikelihoodForSampledAncestors();
    }

}
