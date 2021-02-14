package testdrivers;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.StandardData;
import beast.evolution.datatype.UserDataType;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.prunelikelihood.OrderedDiscreteTraits;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class OrderedDiscreteTraitsTestDriver {

    public static void main(String[] args) {
        // tree
        String treeStr = "((Dog:1.0,Bear:1.0):1.0,Cat:2.0);";
        Tree tree = new TreeParser(treeStr, false, false, true, 0);

        RealParameter liabilities  = new RealParameter(new Double[] {0.0});
        RealParameter thresholds  = new RealParameter(new Double[] {0.0});

        List<Sequence> sequenceList = new ArrayList<>(3);
        Sequence sp1 = new Sequence("Dog", "134");
        Sequence sp2 = new Sequence("Cat", "211");
        Sequence sp3 = new Sequence("Bear", "022");
        sequenceList.add(0, sp1);
        sequenceList.add(1, sp2);
        sequenceList.add(2, sp3);

        UserDataType userDataType1 = new UserDataType();
        userDataType1.initByName("characterName", "ch1","codeMap", "0=0, 1=1, 2=2", "states", 3, "value", "0 no toes, 1 one toe, 2 two toes");
        UserDataType userDataType2 = new UserDataType();
        userDataType2.initByName("characterName", "ch2","codeMap", "0=0, 1=1, 2=2, 3=3", "states", 4, "value", "0 not covered, 1 minor covered, 2 major covered, 3 full covered");
        UserDataType userDataType3 = new UserDataType();
        userDataType3.initByName("characterName", "ch3","codeMap", "0=0, 1=1, 2=2, 3=3, 4=4", "states", 5, "value", "0 complete, 1 a quarter split, 2 half split, 3 three quarters split, 4 full split");

        List<UserDataType> charStateLabels= new ArrayList<>(3);
        charStateLabels.add(0, userDataType1);
        charStateLabels.add(1, userDataType2);
        charStateLabels.add(2, userDataType3);

        StandardData standardData = new StandardData();
        standardData.initByName("charstatelabels", charStateLabels);

        Alignment data = new Alignment();
        data.initByName("sequence", sequenceList, "userDataType", standardData);

        OrderedDiscreteTraits orderedDiscreteTraits = new OrderedDiscreteTraits();
        orderedDiscreteTraits.initByName("liability", liabilities, "threshold", thresholds, "data", data, "tree", tree);

        int[] discreteDataArr = orderedDiscreteTraits.getDiscreteDataArr();
        System.out.println("Array of trait data for all species " + Arrays.toString(discreteDataArr));

        double[] initLiabilities = orderedDiscreteTraits.getLiabilities();
        System.out.println("Initiated liabilities " + Arrays.toString(initLiabilities));

        double[] initThresholds = orderedDiscreteTraits.getThresholds();
        System.out.println("Initiated thresholds " + Arrays.toString(initThresholds));

        orderedDiscreteTraits.calculateLogP();
        double lnLik = orderedDiscreteTraits.getLogP();
        System.out.println("Initial log likelihood = " + lnLik);
    }
}