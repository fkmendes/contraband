package test;

import beast.util.TreeParser;
import contraband.prunelikelihood.MorphologicalData;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Test;
import outercore.parameter.KeyRealParameter;

import java.util.Arrays;
import java.util.List;

public class MorphologicalDataTest {
    private TreeParser tree;
    private String treeStr;
    private String spNames;
    private Integer nTraits;
    private List<Double> contTraitData;
    private List<Double> popTraitData;
    private MorphologicalData morphData = new MorphologicalData();
    private final KeyRealParameter contTrait = new KeyRealParameter();
    private final KeyRealParameter population = new KeyRealParameter();

    @Test
    public void testBMLikelihood() {
        // tree
        treeStr = "";
        spNames = "";
        tree = new TreeParser(treeStr, false, false, true, 0);

        // trait values
        nTraits = 2;
        contTraitData =  Arrays.asList(
        );
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", nTraits);
        morphData.initByName("traits", contTrait, "tree", tree);


        // trait values
        nTraits = 2;
        RealMatrix traitValuesRM = new Array2DRowRealMatrix(new double[][]{
                {-2.62762948691895, -1.56292164859448},
                {-1.50846427625826, -1.59482814741543},
                {-0.226074849617958, -2.11000367246907}
        });

        // shrinkage parameter for the population variance
        double lambda = 0.574732079259225;
    }

}
