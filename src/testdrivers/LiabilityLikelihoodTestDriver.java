package testdrivers;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import contraband.clock.RateCategoryClockModel;
import contraband.math.MatrixUtilsContra;
import contraband.math.NodeMath;
import contraband.prunelikelihood.LiabilityLikelihood;
import outercore.parameter.KeyRealParameter;
import java.util.Arrays;
import java.util.List;

public class LiabilityLikelihoodTestDriver {

    public static void main(String[] args) {
        // tree
        String treeStr = "(((t4:0.01700026351,t10:0.01700026351):0.9034244971,(t2:0.2519264674,t1:0.2519264674):0.6684982932):0.2901851714,((t9:0.4243929597,t5:0.4243929597):0.428956968,((t3:0.1329259849,t6:0.1329259849):0.5897285177,(t7:0.3444839115,t8:0.3444839115):0.3781705911):0.1306954251):0.3572600043):0.0;";
        String spNames = "t4 t10 t2 t1 t9 t5 t3 t6 t7 t8";
        Tree tree = new TreeParser(treeStr, false, false, true, 0);

        // liabilities for binary discrete traits
        // three traits for each of the three species
        // NOTE: the liabilities are in order of tip nodes
        // for example, first five elements are liabilities of node 0 which corresponds to species t1
        // and second five elements are liabilities of node 1 which corresponds to species t10
        // and third five elements are liabilities of node 2 which corresponds to species t2
        RealParameter liabilities = new RealParameter(new Double[] {
                -2.16553733680629, -0.093001131322756, -3.32603366513231, 1.53951108454289, -2.64755502918416,
                3.3252839398938, -4.21928449405138, 0.150249836685309, -4.97807077472762, 1.03044047843148,
                0.345250157267404, -0.749225468028548, -3.41701294175762, 0.994478698194049, -3.44063534931992,
                3.60403255588672, 1.3112428586143, -0.732024175817608, -3.42667016024356, 2.3074920977748,
                4.09620039236067, -5.17365781499908, 0.0830324263435673, -3.92697498145996, 1.23323556462526,
                1.77512813315887, 6.02316315859529, -2.73207211765912, -3.7025682293389, -0.571614944120547,
                2.25231513022396, 3.09693981602772, -1.64956110239375, -2.18351165272075, 2.14200618566428,
                8.90155578270267, 0.0936214074138489, -3.00924775897672, -2.18342919122651, 1.66729459543154,
                6.47320619083183, 3.74483134602231, -3.62473323746156, -5.92016122137718, 3.14670114612779,
                2.99774158156863, 3.65399830020205, -1.69746045226351, -5.82006681351925, 0.84039105047492
        });

        // continuous trait values
        // two traits for each of the three species
        int contTraitNr = 5;
        List<Double> contTraitData =  Arrays.asList(
                -2.28731878531972, 8.16016797286394, -1.60790363996585, -3.50375146318285, 2.4180310750393, -0.811378577256192, 8.18829611082131, -2.45648222243076, -4.4274634946743, 3.03356006715664, 13.6231428892787, -8.90573230225415, 4.21682244085758, 4.60332633627662, -0.716172061308149, 7.16152550409223, -12.1514214682735, 5.53944577979937, 7.14147083928846, -1.89982704153151, -0.519423079492719, -5.77669193046342, -0.312497972990135, -1.59583025877346, 1.3446283692346, -0.984130963068706, -3.0397054776387, 6.09565103724857, -0.840066990821478, 3.15088198649174, 6.92575877256499, -6.37840622885572, -1.76067360986415, -3.77388932553854, -0.453482143262971, 5.0957217584282, -5.25782891756779, 2.13664343352796, -2.64427540812803, -2.19363406695927, 5.26074139751262, 7.6536042206705, 5.13405841572984, -8.23652751442488, -0.431284906380698, 4.83917711885396, 1.77623369001487, 3.99553867000827, -10.5516684942347, -3.257852292301
        );
        KeyRealParameter contTrait = new KeyRealParameter();
        contTrait.initByName("value", contTraitData, "keys", spNames, "minordimension", contTraitNr);

        // model parameters
        NodeMath nodeMath = new NodeMath();
        RealParameter sigmasq = new RealParameter(new Double[]{
                3.92444342735533, 3.02812661102414, 3.06556360045831, 2.81005713192123, 2.30077876608194, 4.64631271741818, 3.15615698487533, 1.50684016606037, 3.3548490285494, 2.35882255622108
        });
        RealParameter correlation = new RealParameter(new Double[]{
                -0.714399955235422, -0.170907328370959, -0.172551347408444, -0.262309098150581, -0.695110504515469, -0.722387873101979, -0.533931801095605, -0.0680750994943082, -0.468054719269276, 0.715655430685729, -0.908337666653097, -0.11559985158965, 0.597849691286683, -0.75620148004964, 0.121895967517048, -0.58693722076714, -0.744936699513346, 0.506615728605539, 0.790090718306601, -0.251074448227882, 0.33023038925603, -0.810318678151816, -0.232060724403709, -0.451232710853219, 0.629280077759176, -0.102967317216098, 0.620128706097603, 0.624779019039124, 0.588684642221779, -0.120336624793708, 0.508950317278504, 0.258442263118923, 0.420364802703261, -0.998750453349203, -0.0493668518029153, -0.55976222967729, -0.24036692455411, 0.225542006548494, -0.296404181513935, -0.777729151304811, -0.512761054560542, 0.336111174896359, -0.164706440642476, 0.576391668058932, -0.794270711485296
        });
        nodeMath.initByName( "nTraits", 10, "nSpecies", 10, "sigmasq", sigmasq, "correlation", correlation, "upperMatrix", true);

        RateCategoryClockModel lsc = new RateCategoryClockModel();
        IntegerParameter colorAssignments = new IntegerParameter(new Integer[] {0});
        RealParameter colorValues = new RealParameter(new Double[] {1.0});
        lsc.initByName("nCat", 1, "rateCatAssign", colorAssignments, "rates", colorValues, "tree", tree);

        LiabilityLikelihood liabilityLikelihood = new LiabilityLikelihood();
        liabilityLikelihood.initByName("binaryLiability", liabilities, "traits", contTrait,
                "tree", tree, "nodeMath", nodeMath, "branchRateModel", lsc);

        double[] combinedDataArr = liabilityLikelihood.getCombinedTraitDataArr();
        System.out.println("Array of combined trait values: " + "\n" + Arrays.toString(combinedDataArr));

        double lik = liabilityLikelihood.calculateLogP();
        System.out.println("Log likelihood"+ lik);
        // expected: -205.515395500167
    }

}
