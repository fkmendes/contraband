package testdrivers;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.StandardData;
import beast.evolution.datatype.UserDataType;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.prunelikelihood.UnorderedDiscreteTraits;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class UnorderedDiscreteTraitsTestDriver {

    public static void main(String[] args) {
        // tree
        String treeStr = "((Dog:1.0,Bear:1.0):1.0,Cat:2.0);";
        Tree tree = new TreeParser(treeStr, false, false, true, 0);

        RealParameter liabilities  = new RealParameter(new Double[] {0.0});

        List<Sequence> sequenceList = new ArrayList<>(3);
        Sequence sp1 = new Sequence("Dog", "101");
        Sequence sp2 = new Sequence("Cat", "243");
        Sequence sp3 = new Sequence("Bear", "032");
        sequenceList.add(0, sp1);
        sequenceList.add(1, sp2);
        sequenceList.add(2, sp3);

        UserDataType userDataType1 = new UserDataType();
        userDataType1.initByName("characterName", "ch1","codeMap", "0=0, 1=1, 2=2", "states", 3, "value", "0 red, 1 blue, 2 yellow");
        UserDataType userDataType2 = new UserDataType();
        userDataType2.initByName("characterName", "ch2","codeMap", "0=0, 1=1, 2=2, 3=3, 4=4", "states", 5, "value", "0 round, 1 triangular, 3 square, 4 rectangular");
        UserDataType userDataType3 = new UserDataType();
        userDataType3.initByName("characterName", "ch3","codeMap", "0=0, 1=1, 2=2, 3=3", "states", 4, "value", "0 two wings, 1 four legs, 2 six paws, 3 two forelimbs");

        List<UserDataType> charStateLabels= new ArrayList<>(3);
        charStateLabels.add(0, userDataType1);
        charStateLabels.add(1, userDataType2);
        charStateLabels.add(2, userDataType3);

        StandardData standardData = new StandardData();
        standardData.initByName("charstatelabels", charStateLabels);

        Alignment data = new Alignment();
        data.initByName("sequence", sequenceList, "userDataType", standardData);

        UnorderedDiscreteTraits unorderedDiscreteTraits = new UnorderedDiscreteTraits();
        unorderedDiscreteTraits.initByName("liability", liabilities, "data", data, "tree", tree);

        int[] discreteDataArr = unorderedDiscreteTraits.getDiscreteDataArr();
        System.out.println("Array of trait data for all species " + Arrays.toString(discreteDataArr));

        double[] initLiabilities = unorderedDiscreteTraits.getLiabilities();
        System.out.println("Initiated liabilities " + Arrays.toString(initLiabilities));

        unorderedDiscreteTraits.calculateLogP();
        double lnLik = unorderedDiscreteTraits.getLogP();
        System.out.println("Initial log likelihood = " + lnLik);
    }
}
