package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.StandardData;
import beast.evolution.datatype.UserDataType;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.prunelikelihood.BinaryDiscreteTraits;
import contraband.prunelikelihood.LiabilityLikelihood;
import contraband.prunelikelihood.OrderedDiscreteTraits;
import contraband.prunelikelihood.UnorderedDiscreteTraits;
import org.omg.PortableInterceptor.SYSTEM_EXCEPTION;
import outercore.parameter.KeyRealParameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LiabilityLikelihoodTestDriver {

    public static void main(String[] args) {
        // tree
        String treeStr = "(((t4:0.01700026351,t10:0.01700026351):0.9034244971,(t2:0.2519264674,t1:0.2519264674):0.6684982932):0.2901851714,((t9:0.4243929597,t5:0.4243929597):0.428956968,((t3:0.1329259849,t6:0.1329259849):0.5897285177,(t7:0.3444839115,t8:0.3444839115):0.3781705911):0.1306954251):0.3572600043):0.0;";
        String spNames = "t4 t10 t2 t1 t9 t5 t3 t6 t7 t8";
        Tree tree = new TreeParser(treeStr, false, false, true, 0);

        List<Sequence> sequenceList = new ArrayList<>(10);
        Sequence sp1 = new Sequence("t4", "001034101");
        Sequence sp2 = new Sequence("t10", "010032243");
        Sequence sp3 = new Sequence("t2", "000123032");
        Sequence sp4 = new Sequence("t1", "111123120");
        Sequence sp5 = new Sequence("t9", "110110110");
        Sequence sp6 = new Sequence("t5", "110201110");
        Sequence sp7 = new Sequence("t3", "101111002");
        Sequence sp8 = new Sequence("t6", "010010202");
        Sequence sp9 = new Sequence("t7", "110001041");
        Sequence sp10 = new Sequence("t8", "001101211");
        sequenceList.add(0, sp1);
        sequenceList.add(1, sp2);
        sequenceList.add(2, sp3);
        sequenceList.add(3, sp4);
        sequenceList.add(4, sp5);
        sequenceList.add(5, sp6);
        sequenceList.add(6, sp7);
        sequenceList.add(7, sp8);
        sequenceList.add(8, sp9);
        sequenceList.add(9, sp10);

        UserDataType userDataType1 = new UserDataType();
        userDataType1.initByName("characterName", "ch1","codeMap", "0=0, 1=1", "states", 2, "value", "0 Present, 1 absent");
        UserDataType userDataType2 = new UserDataType();
        userDataType2.initByName("characterName", "ch2","codeMap", "0=0, 1=1", "states", 2, "value", "0 Round, 1 Square");
        UserDataType userDataType3 = new UserDataType();
        userDataType3.initByName("characterName", "ch3","codeMap", "0=0, 1=1", "states", 2, "value", "0 Enlarged, 1 Reduced");

        UserDataType userDataType4 = new UserDataType();
        userDataType4.initByName("characterName", "ch1","codeMap", "0=0, 1=1, 2=2", "states", 3, "value", "0 no toes, 1 one toe, 2 two toes");
        UserDataType userDataType5 = new UserDataType();
        userDataType5.initByName("characterName", "ch2","codeMap", "0=0, 1=1, 2=2, 3=3", "states", 4, "value", "0 not covered, 1 minor covered, 2 major covered, 3 full covered");
        UserDataType userDataType6 = new UserDataType();
        userDataType6.initByName("characterName", "ch3","codeMap", "0=0, 1=1, 2=2, 3=3, 4=4", "states", 5, "value", "0 complete, 1 a quarter split, 2 half split, 3 three quarters split, 4 full split");

        UserDataType userDataType7 = new UserDataType();
        userDataType7.initByName("characterName", "ch1","codeMap", "0=0, 1=1, 2=2", "states", 3, "value", "0 red, 1 blue, 2 yellow");
        UserDataType userDataType8 = new UserDataType();
        userDataType8.initByName("characterName", "ch2","codeMap", "0=0, 1=1, 2=2, 3=3, 4=4", "states", 5, "value", "0 round, 1 triangular, 3 square, 4 rectangular");
        UserDataType userDataType9 = new UserDataType();
        userDataType9.initByName("characterName", "ch3","codeMap", "0=0, 1=1, 2=2, 3=3", "states", 4, "value", "0 two wings, 1 four legs, 2 six paws, 3 two forelimbs");

        List<UserDataType> charStateLabels= new ArrayList<>(9);
        charStateLabels.add(0, userDataType1);
        charStateLabels.add(1, userDataType2);
        charStateLabels.add(2, userDataType3);
        charStateLabels.add(3, userDataType4);
        charStateLabels.add(4, userDataType5);
        charStateLabels.add(5, userDataType6);
        charStateLabels.add(6, userDataType7);
        charStateLabels.add(7, userDataType8);
        charStateLabels.add(8, userDataType9);
        StandardData standardData = new StandardData();
        standardData.initByName("charstatelabels", charStateLabels);
        Alignment data = new Alignment();
        data.initByName("sequence", sequenceList, "userDataType", standardData);


        FilteredAlignment data1 = new FilteredAlignment();
        data1.initByName("data", data, "filter", "1,2,3");
        FilteredAlignment data2 = new FilteredAlignment();
        data2.initByName("data", data, "filter", "4,5,6");
        FilteredAlignment data3 = new FilteredAlignment();
        data3.initByName("data", data, "filter", "7,8,9");

        // liabilities for discrete traits
        RealParameter binaryLiabilities = new RealParameter(new Double[] {0.0});
        RealParameter orderedLiabilities = new RealParameter(new Double[] {0.0});
        RealParameter thresholds  = new RealParameter(new Double[] {0.0});
        RealParameter unorderedLiabilities  = new RealParameter(new Double[] {0.0});

        // initiate discrete traits
        BinaryDiscreteTraits binaryDiscreteTraits = new BinaryDiscreteTraits();
        binaryDiscreteTraits.initByName("liability", binaryLiabilities, "data", data1, "tree", tree);
        System.out.println("Print binary discrete data " + "\n" + Arrays.toString(binaryDiscreteTraits.getDiscreteDataArr()));
        System.out.println("Print initial binary liabilities " + "\n" + Arrays.toString(binaryDiscreteTraits.getLiabilities()));
        System.out.println("Print likelihood for binary liabilities " + binaryDiscreteTraits.calculateLogP());

        OrderedDiscreteTraits orderedDiscreteTraits = new OrderedDiscreteTraits();
        orderedDiscreteTraits.initByName("liability", orderedLiabilities, "data", data2, "threshold", thresholds, "tree", tree);
        System.out.println("Print ordered discrete data " + "\n" + Arrays.toString(orderedDiscreteTraits.getDiscreteDataArr()));
        System.out.println("Print initial ordered liabilities " + "\n" + Arrays.toString(orderedDiscreteTraits.getLiabilities()));
        System.out.println("Print initial thresholds " + "\n" + Arrays.toString(orderedDiscreteTraits.getThresholds()));
        System.out.println("Print likelihood for ordered liabilities " + orderedDiscreteTraits.calculateLogP());

        UnorderedDiscreteTraits unorderedDiscreteTraits = new UnorderedDiscreteTraits();
        unorderedDiscreteTraits.initByName("liability", unorderedLiabilities, "data", data3, "tree", tree);
        System.out.println("Print unordered discrete data " + "\n" + Arrays.toString(unorderedDiscreteTraits.getDiscreteDataArr()));
        System.out.println("Print initial unordered liabilities " + "\n" + Arrays.toString(unorderedDiscreteTraits.getLiabilities()));
        System.out.println("Print likelihood for unordered liabilities " + unorderedDiscreteTraits.calculateLogP());

        // continuous trait values
        int contTraitNr = 5;
        List<Double> contTraitData =  Arrays.asList(
                -2.28731878531972, 8.16016797286394, -1.60790363996585, -3.50375146318285, 2.4180310750393, -0.811378577256192, 8.18829611082131, -2.45648222243076, -4.4274634946743, 3.03356006715664, 13.6231428892787, -8.90573230225415, 4.21682244085758, 4.60332633627662, -0.716172061308149, 7.16152550409223, -12.1514214682735, 5.53944577979937, 7.14147083928846, -1.89982704153151, -0.519423079492719, -5.77669193046342, -0.312497972990135, -1.59583025877346, 1.3446283692346, -0.984130963068706, -3.0397054776387, 6.09565103724857, -0.840066990821478, 3.15088198649174, 6.92575877256499, -6.37840622885572, -1.76067360986415, -3.77388932553854, -0.453482143262971, 5.0957217584282, -5.25782891756779, 2.13664343352796, -2.64427540812803, -2.19363406695927, 5.26074139751262, 7.6536042206705, 5.13405841572984, -8.23652751442488, -0.431284906380698, 4.83917711885396, 1.77623369001487, 3.99553867000827, -10.5516684942347, -3.257852292301
        );
        KeyRealParameter contTrait = new KeyRealParameter();
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", contTraitNr);

        // model parameters
        NodeMath nodeMath = new NodeMath();
        RealParameter sigmasq = new RealParameter(new Double[]{
                2.35882255622108, 1.97319269488477, 2.54621393782772, 1.99810386929218, 2.18438130143428, 2.25350823993028, 1.63884546108462, 3.49500956814383, 2.84627731546465, 1.93200614977075, 3.95960580579731, 3.08928374542405, 2.48799845800124, 3.55564932692539, 3.5375699840826, 3.47805880846592, 3.34207568397888, 3.20969300681193, 2.66825976585302
        });
        RealParameter correlation = new RealParameter(new Double[]{
                -0.24036692455411, 0.225542006548494, -0.296404181513935, -0.777729151304811, -0.512761054560542, 0.336111174896359, -0.164706440642476, 0.576391668058932, -0.794270711485296, -0.130214517004788, 0.969913959968835, 0.786102228797972, 0.772938121575862, -0.649894699454308, -0.738608616869897, 0.306203850079328, -0.312967055477202, 0.313516255933791, -0.359253515023738, -0.624617761466652, 0.564588602632284, -0.81281002657488, -0.06644191686064, 0.0230109198018909, 0.19997791852802, -0.334352919366211, -0.0227739326655865, 0.908947654999793, -0.0341952056623995, 0.780700444243848, 0.828876373823732, 0.217469964642078, -0.178620446939021, -0.705810618121177, 0.870599606540054, -0.397542200051248, -0.878558856900781, 0.895453880075365, 0.441192546859384, -0.715411408804357, 0.098569312132895, 0.908182477112859, 0.170966706238687, -0.190979436505586, 0.29578695865348, -0.360358765814453, -0.384559978265315, -0.560464737471193, -0.261022268328816, 0.96843840694055, -0.691595398355275, -0.817912000231445, -0.716186184436083, 0.380014203023165, 0.238512966781855, 0.782788234297186, 0.345998185221106, 0.474155475851148, 0.0422714515589178, 0.31967689935118, 0.643610920291394, 0.572563103400171, 0.959643834736198, -0.121136927511543, -0.376595595851541, -0.181050094775856, -0.979065776336938, -0.63230095198378, 0.685458637773991, -0.537676435895264, -0.521800088696182, -0.846617669332772, -0.508552643936127, 0.46427041105926, 0.694906330201775, -0.00494546582922339, -0.224181940313429, -0.507102011702955, -0.777807077392936, -0.220011129509658, 0.143870627973229, -0.566214474383742, -0.110463995952159, -0.564018662553281, 0.00459912652149796, -0.292190856300294, 0.29997031763196, -0.250572086777538, -0.289109238423407, 0.0673758909106255, 0.480668720789254, -0.557794124353677, -0.1745077627711, -0.468626626301557, 0.259946106933057, -0.632343018427491, 0.727288222871721, 0.493136008270085, 0.336569299455732, 0.236035746522248, -0.255523879546672, 0.0596713717095554, 0.749364685732871, 0.163500199560076, 0.679535529576242, -0.375103670172393, 0.416580644436181, -0.469964387826622, 0.188686388079077, -0.0374203990213573, -0.46993453707546, 0.129180869553238, 0.826376446057111, 0.803748778998852, -0.451666756998748, -0.357034487184137, 0.971281768754125, 0.239986620377749, 0.874628178309649, -0.0669345953501761, -0.186334813479334, 0.318460648413748, -0.695306766312569, 0.145734116435051, -0.522547946311533, 0.924717872869223, 0.202731451950967, 0.0300594544969499, -0.194853315595537, 0.760493082460016, -0.271816270425916, -0.423521438613534, -0.658709529787302, -0.655656507238746, -0.0359147889539599, -0.494070142973214, -0.5674904207699, 0.348752776160836, -0.904672745149583, 0.401706174947321, -0.296222723089159, -0.182112004142255, 0.641902647912502, 0.83771469630301, -0.434943339787424, 0.922209587413818, 0.45678885653615, 0.372750164009631, -0.894312114454806, -0.209559730719775, -0.0443092403002083, 0.120506527367979, 0.396523189730942, 0.831367076840252, 0.236702454742044, -0.143156982492656, 0.0841607344336808, -0.883043022826314, -0.478286285884678, -0.205696093384176, -0.604510526638478, 0.663855125661939, -0.694225554354489, 0.606837084051222, 0.0936523131094873, 0.324635284021497, -0.656603012233973, 0.266110719647259, -0.376260506454855, 0.44910869281739, -0.202120350673795
        });
        nodeMath.initByName( "binaryDiscreteTraits", binaryDiscreteTraits,
                "orderedDiscreteTraits", orderedDiscreteTraits,
                "traits", contTrait,
                "unorderedDiscreteTraits", unorderedDiscreteTraits, "sigmasq", sigmasq, "correlation", correlation, "upperMatrix", true);

        RateCategoryClockModel lsc = new RateCategoryClockModel();
        IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
        RealParameter colorValues = new RealParameter(new Double[] {1.0});
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        LiabilityLikelihood liabilityLikelihood = new LiabilityLikelihood();
        liabilityLikelihood.initByName("binaryDiscreteTraits", binaryDiscreteTraits,
                "orderedDiscreteTraits", orderedDiscreteTraits,
                "traits", contTrait,
                "unorderedDiscreteTraits", unorderedDiscreteTraits,
                "tree", tree, "nodeMath", nodeMath, "branchRateModel", lsc);

        double[] combinedDataArr = liabilityLikelihood.getCombinedTraitDataArr();
        System.out.println("Array of combined trait values: " + "\n" + Arrays.toString(combinedDataArr));

        double lik = liabilityLikelihood.calculateLogP();
        System.out.println("Log likelihood"+ lik);//
        // expected:
    }

}
