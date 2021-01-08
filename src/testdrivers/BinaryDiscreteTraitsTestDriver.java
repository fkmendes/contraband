package testdrivers;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.StandardData;
import beast.evolution.datatype.UserDataType;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.prunelikelihood.BinaryDiscreteTraits;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BinaryDiscreteTraitsTestDriver {

    public static void main(String[] args) {
        // tree
        String treeStr = "((Dog:1.0,Bear:1.0):1.0,Cat:2.0);";
        Tree tree = new TreeParser(treeStr, false, false, true, 0);

        RealParameter liabilities  = new RealParameter(new Double[] {0.0});

        List<Sequence> sequenceList = new ArrayList<>(3);
        Sequence sp1 = new Sequence("Dog", "1100");
        Sequence sp2 = new Sequence("Cat", "0111");
        Sequence sp3 = new Sequence("Bear", "0010");
        sequenceList.add(0, sp1);
        sequenceList.add(1, sp2);
        sequenceList.add(2, sp3);

        UserDataType userDataType1 = new UserDataType();
        userDataType1.initByName("characterName", "ch1","codeMap", "0=0, 1=1", "states", 2, "value", "0 Present, 1 absent");
        UserDataType userDataType2 = new UserDataType();
        userDataType2.initByName("characterName", "ch2","codeMap", "0=0, 1=1", "states", 2, "value", "0 Round, 1 Square");
        UserDataType userDataType3 = new UserDataType();
        userDataType3.initByName("characterName", "ch3","codeMap", "0=0, 1=1", "states", 2, "value", "0 Enlarged, 1 Reduced");

        List<UserDataType> charStateLabels= new ArrayList<>(3);
        charStateLabels.add(0, userDataType1);
        charStateLabels.add(1, userDataType2);
        charStateLabels.add(2, userDataType3);

        StandardData standardData = new StandardData();
        standardData.initByName("charstatelabels", charStateLabels);

        Alignment data = new Alignment();
        data.initByName("sequence", sequenceList, "userDataType", standardData);

        BinaryDiscreteTraits discreteTraits = new BinaryDiscreteTraits();
        discreteTraits.initByName("liability", liabilities, "data", data, "tree", tree);

        int[] sp1Data = discreteTraits.getAllDataForSpecies(tree.getNode(0).getID());
        System.out.println("Trait for " + tree.getNode(0).getID() + " is " + Arrays.toString(sp1Data));

        int[] sp2Data = discreteTraits.getAllDataForSpecies(tree.getNode(1).getID());
        System.out.println("Trait for " + tree.getNode(1).getID() + " is " + Arrays.toString(sp2Data));

        int[] sp3Data = discreteTraits.getAllDataForSpecies(tree.getNode(2).getID());
        System.out.println("Trait for " + tree.getNode(2).getID() + " is " + Arrays.toString(sp3Data));

        int[] discreteDataArr = discreteTraits.getDiscreteDataArr();
        System.out.println("Array of trait data for all species " + Arrays.toString(discreteDataArr));

        double[] initLiabilities = discreteTraits.getLiabilities();
        System.out.println("Initiated liabilities " + Arrays.toString(initLiabilities));

        discreteTraits.calculateLogP();

        double lnLik = discreteTraits.getLogP();
        System.out.println("Initial log likelihood = " + lnLik);
    }
}
