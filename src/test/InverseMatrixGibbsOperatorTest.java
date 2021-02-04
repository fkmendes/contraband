package test;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeTrait;
import beast.math.distributions.WishartDistribution;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.NodeMath;
import contraband.prunelikelihood.InverseMatrixGibbsOperator;
import contraband.prunelikelihood.LiabilityLikelihood;
import org.junit.Assert;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

public class InverseMatrixGibbsOperatorTest {

    private final static double EPSILON = 1e-6;
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> data;
    private RealParameter inverseMatrix;
    private RealParameter sigmasq;
    private RealParameter correlation;
    private RealParameter rootValues;
    private final NodeMath nodeMath = new NodeMath();
    private final KeyRealParameter traitValues = new KeyRealParameter();
    private final IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
    private final RealParameter colorValues = new RealParameter(new Double[] {1.0});
    private final RateCategoryClockModel lsc = new RateCategoryClockModel();

    @Test
    public void testOperator () {

        // tree with sampled ancestors
        treeStr = "((t6_1:3.271445528,(t8_2:0.3283657173,t8_1:0):2.943079811):2.339286241,((t3_1:0.9555596829,(t2_2:0.5475127318,t2_1:0):0.408046951):2.320853918,t14_1:0):2.334318168):2.275567552;";
        tree = new TreeParser(treeStr, false, false, true, 0);
        spNames = "t8_2 t6_1 t3_1 t2_2 t8_1 t2_1 t14_1";

        // trait values
        nTraits = 6;
        data = Arrays.asList(
                -1.73589276210399, 2.73121128831183, -3.92957539479266, -5.92955994852884, 4.42485042137109, -3.3483724622711, -1.98359780509439, 4.67446542773907, -3.19870142692452, -2.48575462626222, -2.15377773070886, -1.78421025088891, -2.48459828129598, 0.989252696366672, 3.91272385441485, -4.28129680647532, -0.790392197650521, 2.83899780891758, -3.09710366090904, 1.8305494605968, 0.855740412235961, -4.21385872055373, 0.397647886387143, 1.00031968784065, -1.87658975668587, 2.04949154345486, -3.52443574180981, -6.79449126307887, 5.31605700662329, -3.00850971793145, -3.0530711014376, 2.40586268291689, 1.34517681660018, -3.14042322388396, -1.59735551807305, 1.57493293491501, -2.85676462516979, 2.02192686323489, 1.09419702953616, -3.76587499531895, -0.602026063527102, 1.27919157725599
        );
        traitValues.initByName("value", data, "keys", spNames, "minordimension", nTraits);

        // trait evolutionary rate matrix -> inverse matrix
        sigmasq = new RealParameter(new Double[]{1.0});
        inverseMatrix = new RealParameter(new Double[] {20.0860540472875, -0.699178622772761, -11.1260163333518, -7.38252183916494, -3.16775577164432, 15.4353824245299, -0.699178622772761, 39.9308170368297, 22.6721425407894, 7.66130953107541, 21.6872726722536, 0.230109078471028, -11.1260163333518, 22.6721425407894, 22.7027211555621, 7.71928485251703, 12.0533069136996, -15.0184390313218, -7.38252183916494, 7.66130953107541, 7.71928485251703, 11.4886127503156, 9.41044618131386, -1.60210698031551, -3.16775577164432, 21.6872726722536, 12.0533069136996, 9.41044618131386, 15.6478177926314, 2.43677913653119, 15.4353824245299, 0.230109078471028, -15.0184390313218, -1.60210698031551, 2.43677913653119, 25.9088965100293});
        nodeMath.initByName("traits", traitValues, "sigmasq", sigmasq, "inverseMatrix", inverseMatrix);

        // branch rate model
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        // likelihood
        LiabilityLikelihood likelihood = new LiabilityLikelihood();
        likelihood.initByName("nodeMath", nodeMath, "tree", tree, "traits", traitValues, "branchRateModel", lsc);

        RealParameter scaleMatrix = new RealParameter(new Double[] {6.0, 1.74420791678131, 0.393594525754452, -3.25641066860408, -1.13246589247138, 0.178297321312129, 1.74420791678131, 6.0, 4.30627778358757, 0.638826004229484, 0.421550226397813, -1.1075679268688, 0.393594525754452, 4.30627778358757, 6.0, 1.95473941415548, -1.3254970824346, -3.68255945853889, -3.25641066860408, 0.638826004229484, 1.95473941415548, 6.0, 1.28835505899042, 0.316137732937932, -1.13246589247138, 0.421550226397813, -1.3254970824346, 1.28835505899042, 6.0, 0.69566449150443, 0.178297321312129, -1.1075679268688, -3.68255945853889, 0.316137732937932, 0.69566449150443, 6.0});
        WishartDistribution distribution = new WishartDistribution();
        distribution.initByName("df", 6, "scaleMatrix", scaleMatrix);

        RealParameter values = new RealParameter(new Double[]{1.0});
        InverseMatrixGibbsOperator operator = new InverseMatrixGibbsOperator();
        operator.initByName("likelihood", likelihood, "nodeMath", nodeMath,
                "dimension", 4,
                "distr", distribution,
                "values", values);

        operator.proposal();

    }
}
