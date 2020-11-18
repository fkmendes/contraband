package testdrivers;

import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import contraband.prunelikelihood.OUPruneLikelihood;
import contraband.prunelikelihood.OUPruneUtils;
import contraband.utils.GeneralUtils;
import contraband.valuewrappers.OneValueContTraits;
import org.apache.commons.math3.linear.*;
import java.util.ArrayList;
import java.util.List;
import beast.evolution.tree.Node;

public class OUPruneLikelihoodTestDriver3 {
    public static void main(String[] args) {
        String treeStr = "((A:3.0058179,B:3.0058179):4.350951,C:7.3567689);";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        int nTraits = 2;

        RealParameter data = new RealParameter(new Double[] {
                0.131394584822684, -0.65643500381027,
                -0.17399894787897, -0.19269144948362,
                -0.373259036803319, 1.27056761078824
        });
        String spNames = "A,B,C";
        OneValueContTraits traitValues = new OneValueContTraits();
        traitValues.initByName("nTraits", nTraits, "spNames", spNames, "traitValues", data);

        // OU model parameters
        RealParameter rootValues = new RealParameter(new Double[] {0.2, 1.3});

        RealParameter sigma = new RealParameter(new Double[] {1.0, 0.6, 1.0});

        RealParameter sigmae = new RealParameter(new Double[] {1.0, 0.3, 1.0});

        RealParameter alpha = new RealParameter(new Double[] {2.0, 0.0, 0.0, 2.0});

        RealParameter theta = new RealParameter(new Double[] {0.5, 0.5});

        OUPruneLikelihood pcm = new OUPruneLikelihood();
        int nSpecies = tree.getLeafNodeCount();
        int nodeCount = tree.getNodeCount();

        List<RealVector> traitsValuesList = new ArrayList<>(nSpecies);
        pcm.populateTraitsValuesList(traitValues, tree, traitsValuesList);

        RealMatrix sigmaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        pcm.populateSigmaMatrix(sigmaRM, sigma.getDoubleValues());

        RealMatrix alphaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        pcm.populateAlphaMatrix(alphaRM, alpha.getDoubleValues());

        RealVector thetaVec = new ArrayRealVector(new double[nTraits]);
        pcm.populateRealVector(thetaVec, theta.getDoubleValues());

        RealMatrix sigmaERM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        pcm.populateSigmaMatrix(sigmaERM, sigmae.getDoubleValues());

        List<RealMatrix> lMatList = new ArrayList<>(nodeCount);
        List<RealVector> mVecList = new ArrayList<>(nodeCount);
        double [] rArr = new double[nodeCount];
        RealVector rootValuesVec = new ArrayRealVector(rootValues.getDoubleValues());
        double []vcvMatDetArr = new double[nodeCount];
        double [] negativeTwoAplusLDetArr = new double[nodeCount];

        // initiate
        RealMatrix iniLMat = new Array2DRowRealMatrix(new double [nTraits][nTraits]);
        RealVector iniMVec = new ArrayRealVector(new double [nTraits]);
        for (int i= 0; i < tree.getNodeCount(); i ++) {
            lMatList.add(i, iniLMat);
            mVecList.add(i, iniMVec);
            rArr[i] = 0;

        }

        RealMatrix identity = MatrixUtils.createRealIdentityMatrix(nTraits);

        EigenDecomposition decompositionH = new EigenDecomposition(alphaRM);
        RealMatrix pMat = decompositionH.getV();
        RealMatrix inverseP = new LUDecomposition(pMat).getSolver().getInverse();
        //System.out.println("Display H :");
        //GeneralUtils.displayRealMatrix(alphaRM);

        // TEST1: pruning process
        // at Node A
        System.out.println("At Node A :");
        int childIdx = 0;
        Node child = tree.getNode(childIdx);
        RealMatrix thisNodeLMat = iniLMat;
        RealVector thisNodeMVec = iniMVec;
        double thisNodeR = 0;

        RealMatrix phiRM = OUPruneUtils.getPhiRM(child, alphaRM);
        System.out.println("Display phiRM :");
        GeneralUtils.displayRealMatrix(phiRM);
        //0.002450077 0.000000000
        //0.000000000 0.002450077

        RealVector omegaVec = OUPruneUtils.getOmegaVec(thetaVec, phiRM, identity);
        System.out.println("Display megaVec  :");
        GeneralUtils.displayRealVector(omegaVec);
        // 0.498775
        // 0.498775

        RealMatrix varianceRM = OUPruneUtils.getOUVarianceRM(child, sigmaRM, sigmaERM, pMat, inverseP, decompositionH, nTraits);
        System.out.println("Display varianceRM :");
        GeneralUtils.displayRealMatrix(varianceRM );
        // 1.4299980 0.4499991
        // 0.4499991 1.2499985

        RealMatrix invVCVMat = OUPruneUtils.getInverseVarianceRMForOU(child, vcvMatDetArr, sigmaRM, sigmaERM, pMat, inverseP, decompositionH, nTraits);
        double varianceRMDet = vcvMatDetArr[childIdx];

        RealMatrix aMat = OUPruneUtils.getAMatForOU(invVCVMat);
        RealMatrix eMat = OUPruneUtils.getEMatForOU(phiRM, invVCVMat);
        RealMatrix cMat = OUPruneUtils.getCMatForOU(phiRM, eMat);
        double f = OUPruneUtils.getFforOU(omegaVec, invVCVMat, varianceRMDet, nTraits);
        RealVector bVec = OUPruneUtils.getBVecForOU(invVCVMat, omegaVec);
        RealVector dVec = OUPruneUtils.getDVecForOU(eMat, omegaVec);

        System.out.println("Display A :");
        GeneralUtils.displayRealMatrix(aMat);
        //-0.39432226013635 0.14195590002967257
        //0.1419559000296726 -0.4511045633385124
        System.out.println("Display E :");
        GeneralUtils.displayRealMatrix(eMat);
        //0.0019322398512067107 -6.956057896817779E-4
        //-6.95605789681778E-4 0.002210481888703103
        System.out.println("Display C :");
        GeneralUtils.displayRealMatrix(cMat);
        //-2.3670682713298583E-6 8.521438956353147E-7
        //8.521438956353148E-7 -2.7079254885622673E-6
        System.out.println("Display b:");
        GeneralUtils.displayRealVector(bVec);
        //0.2517480430759149
        //0.30839122525932916
        System.out.println("Display d:");
        GeneralUtils.displayRealVector(dVec);
        //-6.168021063867717E-4
        //-7.555822679177356E-4
        System.out.println("Display f = " + f);
        //-2.2078597649235796
        RealVector traitsVec = traitsValuesList.get(childIdx);
        // 0.1313945848; -0.1926914495
        // set the L matrix
        RealMatrix lMat = OUPruneUtils.getLMatForOULeaf(cMat);
        System.out.println("Display L :");
        GeneralUtils.displayRealMatrix(lMat);
        //-2.3670682713298583E-6 8.521438956353147E-7
        //8.521438956353147E-7 -2.7079254885622673E-6
        thisNodeLMat = thisNodeLMat.add(OUPruneUtils.getLMatForOULeaf(cMat));

        // set r value
        double r = OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f);
        System.out.println("Display r = " + r);
        //-2.2649513417863973
        thisNodeR += OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f);

        // set m vector
        RealVector mVec = OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec);
        System.out.println("Display mVec:");
        GeneralUtils.displayRealVector(mVec);
        //-2.2887896547664136E-4
        //-0.0012729220610447189
        thisNodeMVec = thisNodeMVec.add(OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec));

        // at Node B
        System.out.println("At Node B :");
        childIdx = 1;
        child = tree.getNode(childIdx);
        phiRM = OUPruneUtils.getPhiRM(child, alphaRM);
        System.out.println("Display phiRM :");
        GeneralUtils.displayRealMatrix(phiRM);
        //0.002450077064554477 0.0
        //0.0 0.002450077064554477
        omegaVec = OUPruneUtils.getOmegaVec(thetaVec, phiRM, identity);
        System.out.println("Display megaVec  :");
        GeneralUtils.displayRealVector(omegaVec);
        //0.49877496146772277
        //0.49877496146772277
        varianceRM = OUPruneUtils.getOUVarianceRM(child, sigmaRM, sigmaERM, pMat, inverseP, decompositionH, nTraits);
        System.out.println("Display varianceRM :");
        GeneralUtils.displayRealMatrix(varianceRM );
        //1.4299979590216085 0.44999909956835665
        //0.44999909956835665 1.2499984992805944
        invVCVMat = OUPruneUtils.getInverseVarianceRMForOU(child, vcvMatDetArr, sigmaRM, sigmaERM, pMat, inverseP, decompositionH, nTraits);
        varianceRMDet = vcvMatDetArr[childIdx];

        aMat = OUPruneUtils.getAMatForOU(invVCVMat);
        eMat = OUPruneUtils.getEMatForOU(phiRM, invVCVMat);
        cMat = OUPruneUtils.getCMatForOU(phiRM, eMat);
        f = OUPruneUtils.getFforOU(omegaVec, invVCVMat, varianceRMDet, nTraits);
        bVec = OUPruneUtils.getBVecForOU(invVCVMat, omegaVec);
        dVec = OUPruneUtils.getDVecForOU(eMat, omegaVec);

        System.out.println("Display A :");
        GeneralUtils.displayRealMatrix(aMat);
        //-0.39432226013635 0.14195590002967257
        //0.1419559000296726 -0.4511045633385124
        System.out.println("Display E :");
        GeneralUtils.displayRealMatrix(eMat);
        //0.0019322398512067107 -6.956057896817779E-4
        //-6.95605789681778E-4 0.002210481888703103
        System.out.println("Display C :");
        GeneralUtils.displayRealMatrix(cMat);
        //-2.3670682713298583E-6 8.521438956353147E-7
        //8.521438956353148E-7 -2.7079254885622673E-6
        System.out.println("Display b:");
        GeneralUtils.displayRealVector(bVec);
        //0.2517480430759149
        //0.30839122525932916
        System.out.println("Display d:");
        GeneralUtils.displayRealVector(dVec);
        //-6.168021063867717E-4
        //-7.555822679177356E-4
        System.out.println("Display f = " + f);
        //-2.2078597649235796
        traitsVec = traitsValuesList.get(childIdx);

        // set the L matrix
        lMat = OUPruneUtils.getLMatForOULeaf(cMat);
        System.out.println("Display L :");
        GeneralUtils.displayRealMatrix(lMat);
        //-2.3670682713298583E-6 8.521438956353147E-7
        //8.521438956353147E-7 -2.7079254885622673E-6
        thisNodeLMat = thisNodeLMat.add(OUPruneUtils.getLMatForOULeaf(cMat));

        // set r value
        r = OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f);
        System.out.println("Display r = " + r);
        //-2.651426768086688
        thisNodeR += OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f);

        // set m vector
        mVec = OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec);
        System.out.println("Display mVec:");
        GeneralUtils.displayRealVector(mVec);
        //-0.0016255508334245719
        //-0.0011240446193660335
        thisNodeMVec = thisNodeMVec.add(OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec));

        // set Values at Node D
        lMatList.set(3, thisNodeLMat);
        mVecList.set(3, thisNodeMVec);
        rArr[3] = thisNodeR;

        // at Node C
        System.out.println("At Node C :");
        thisNodeLMat = iniLMat;
        thisNodeMVec = iniMVec;
        thisNodeR = 0;
        childIdx = 2;
        child = tree.getNode(childIdx);
        phiRM = OUPruneUtils.getPhiRM(child, alphaRM);
        System.out.println("Display phiRM :");
        GeneralUtils.displayRealMatrix(phiRM);
        //4.073725149339952E-7 0.0
        //0.0 4.073725149339952E-7
        omegaVec = OUPruneUtils.getOmegaVec(thetaVec, phiRM, identity);
        System.out.println("Display megaVec  :");
        GeneralUtils.displayRealVector(omegaVec);
        //0.49999979631374253
        //0.49999979631374253
        varianceRM = OUPruneUtils.getOUVarianceRM(child, sigmaRM, sigmaERM, pMat, inverseP, decompositionH, nTraits);
        System.out.println("Display varianceRM :");
        GeneralUtils.displayRealMatrix(varianceRM );
        //1.4299999999999435 0.4499999999999751
        //0.4499999999999751 1.2499999999999585
        invVCVMat = OUPruneUtils.getInverseVarianceRMForOU(child, vcvMatDetArr, sigmaRM, sigmaERM, pMat, inverseP, decompositionH, nTraits);
        varianceRMDet = vcvMatDetArr[childIdx];

        aMat = OUPruneUtils.getAMatForOU(invVCVMat);
        eMat = OUPruneUtils.getEMatForOU(phiRM, invVCVMat);
        cMat = OUPruneUtils.getCMatForOU(phiRM, eMat);
        f = OUPruneUtils.getFforOU(omegaVec, invVCVMat, varianceRMDet, nTraits);
        bVec = OUPruneUtils.getBVecForOU(invVCVMat, omegaVec);
        dVec = OUPruneUtils.getDVecForOU(eMat, omegaVec);

        System.out.println("Display A :");
        GeneralUtils.displayRealMatrix(aMat);
        //-0.39432176656152784 0.1419558359621469
        //0.14195583596214686 -0.45110410094638503
        System.out.println("Display E :");
        GeneralUtils.displayRealMatrix(eMat);
        //3.2127169947477074E-7 -1.1565781181091492E-7
        //-1.156578118109149E-7 3.675348241991354E-7
        System.out.println("Display C :");
        GeneralUtils.displayRealMatrix(cMat);
        //-6.543863009607804E-14 2.3557906834587576E-14
        //2.355790683458757E-14 -7.48617928299128E-14
        System.out.println("Display b:");
        GeneralUtils.displayRealVector(bVec);
        //0.25236582779243716
        //0.3091481390457319
        System.out.println("Display d:");
        GeneralUtils.displayRealVector(dVec);
        //-1.0280690195120465E-7
        //-1.2593845489022425E-7
        System.out.println("Display f = " + f);
        //-2.2085477045971365
        traitsVec = traitsValuesList.get(childIdx);

        // set the L matrix
        lMat = OUPruneUtils.getLMatForOULeaf(cMat);
        System.out.println("Display L :");
        GeneralUtils.displayRealMatrix(lMat);
        //-6.543863009607804E-14 2.3557906834587572E-14
        //2.3557906834587572E-14 -7.48617928299128E-14
        thisNodeLMat = thisNodeLMat.add(OUPruneUtils.getLMatForOULeaf(cMat));

        // set r value
        r = OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f);
        System.out.println("Display r = " + r);
        //-2.662606604899639
        thisNodeR += OUPruneUtils.getRForOULeaf(aMat, traitsVec, bVec, f);

        // set m vector
        mVec = OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec);
        System.out.println("Display mVec:");
        GeneralUtils.displayRealVector(mVec);
        //-3.056589092646934E-7
        //3.611637261430301E-7
        thisNodeMVec = thisNodeMVec.add(OUPruneUtils.getMVecForOULeaf(eMat, traitsVec, dVec));

        // at Node D
        System.out.println("At Node D :");
        childIdx = 3;
        child = tree.getNode(childIdx);
        phiRM = OUPruneUtils.getPhiRM(child, alphaRM);
        System.out.println("Display phiRM :");
        GeneralUtils.displayRealMatrix(phiRM);
        //1.6626926590493734E-4 0.0
        //0.0 1.6626926590493734E-4
        omegaVec = OUPruneUtils.getOmegaVec(thetaVec, phiRM, identity);
        System.out.println("Display megaVec  :");
        GeneralUtils.displayRealVector(omegaVec);
        //0.49991686536704755
        //0.49991686536704755
        varianceRM = OUPruneUtils.getOUVarianceRM(child, sigmaRM, sigmaERM, pMat, inverseP, decompositionH, nTraits);
        System.out.println("Display varianceRM :");
        GeneralUtils.displayRealMatrix(varianceRM );
        //0.3399999906005406 0.14999999585317966
        //0.14999999585317966 0.2499999930886328
        invVCVMat = OUPruneUtils.getInverseVarianceRMForOU(child, vcvMatDetArr, sigmaRM, sigmaERM, pMat, inverseP, decompositionH, nTraits);
        varianceRMDet = vcvMatDetArr[childIdx];

        aMat = OUPruneUtils.getAMatForOU(invVCVMat);
        eMat = OUPruneUtils.getEMatForOU(phiRM, invVCVMat);
        cMat = OUPruneUtils.getCMatForOU(phiRM, eMat);
        f = OUPruneUtils.getFforOU(omegaVec, invVCVMat, varianceRMDet, nTraits);
        bVec = OUPruneUtils.getBVecForOU(invVCVMat, omegaVec);
        dVec = OUPruneUtils.getDVecForOU(eMat, omegaVec);

        System.out.println("Display A :");
        GeneralUtils.displayRealMatrix(aMat);
        //-2.000000055290939 1.2000000331745633
        //1.2000000331745633 -2.720000075195677
        System.out.println("Display E :");
        GeneralUtils.displayRealMatrix(eMat);
        //6.65077082006117E-4 -3.990462492036702E-4
        //-3.990462492036702E-4 9.045048315283192E-4
        System.out.println("Display C :");
        GeneralUtils.displayRealMatrix(cMat);
        //-5.529093909767744E-8 3.317456345860646E-8
        //3.317456345860646E-8 -7.519567717284132E-8
        System.out.println("Display b:");
        GeneralUtils.displayRealVector(bVec);
        //0.7998670066999745
        //1.5197473127299514
        System.out.println("Display d:");
        GeneralUtils.displayRealVector(dVec);
        //-1.3299330002558432E-4
        //-2.526872700486103E-4
        System.out.println("Display f = " + f);
        //-1.0313898373589483

        // (aMat + lMat).inverse
        RealMatrix aPlusLInv = OUPruneUtils.getInvAPlusLRM(child, negativeTwoAplusLDetArr, aMat, lMatList.get(childIdx));
        System.out.println("Display aPlusLInv = " + aPlusLInv);
        //-0.6799980001,-0.2999989464
        // -0.2999989464,-0.4999987174

        // determinant of -2 * (aMat + lMat)
        double logDetVNode = negativeTwoAplusLDetArr[childIdx];
        System.out.println("Display logDetVNode = " + logDetVNode);
        // 2.7725936820899033

        // set r value
        RealVector a = mVecList.get(childIdx);
        double b = rArr[childIdx];
        r = OUPruneUtils.getRForOUIntNode(bVec, mVecList.get(childIdx), aPlusLInv, f, rArr[childIdx], nTraits, logDetVNode);
        System.out.println("Display r = " + r);
        //-4.91850561594557
        thisNodeR += OUPruneUtils.getRForOUIntNode(bVec, mVecList.get(childIdx), aPlusLInv, f, rArr[childIdx], nTraits, logDetVNode);

        RealMatrix eAPlusLInv = eMat.multiply(aPlusLInv);
        // set m vector
        mVec = OUPruneUtils.getMVecForOUIntNode(eAPlusLInv, bVec, mVecList.get(childIdx), dVec);
        System.out.println("Display mVec:");
        GeneralUtils.displayRealVector(mVec);
        //-3.0883730523252454E-7
        //-3.9915788595441294E-7
        thisNodeMVec = thisNodeMVec.add(OUPruneUtils.getMVecForOUIntNode(eAPlusLInv, bVec, mVecList.get(childIdx), dVec));

        // set L matrix
        lMat = OUPruneUtils.getLMatForOUIntNode(cMat, eMat, eAPlusLInv);
        System.out.println("Display L :");
        GeneralUtils.displayRealMatrix(lMat);
        //-1.3087709637836688E-13 4.711579240003084E-14
        //4.7115792393413396E-14 -1.4972343217692136E-13
        thisNodeLMat = thisNodeLMat.add(OUPruneUtils.getLMatForOUIntNode(cMat, eMat, eAPlusLInv));

        lMatList.set(4, thisNodeLMat);
        mVecList.set(4, thisNodeMVec);
        rArr[4] = thisNodeR;

        // at Node E
        System.out.println("At Node E:");
        System.out.println("Display LE:");
        GeneralUtils.displayRealMatrix(lMatList.get(4));
        //-1.9631572647444492E-13 7.067369923461842E-14
        //7.067369922800097E-14 -2.2458522500683416E-13
        System.out.println("Display mE:");
        GeneralUtils.displayRealVector(mVecList.get(4));
        //-6.144962144972179E-7
        //-3.799415981138285E-8
        System.out.println("Display rE = " + rArr[4]);
        //-7.581112220845209

        // calculate likelihood
        RealMatrix lE = lMatList.get(tree.getRoot().getNr());
        RealVector mE = mVecList.get(tree.getRoot().getNr());
        double rE = rArr[tree.getRoot().getNr()];
        double loglik1 = lE.preMultiply(rootValuesVec).dotProduct(rootValuesVec) + rootValuesVec.dotProduct(mE) + rE;
        System.out.println("Log likelihood1 = " + loglik1);
        // expected: -7.58111239313721

        // TEST2: pruneOUPCM()
        OUPruneUtils.pruneOUPCM(tree.getRoot(), nTraits, traitsValuesList, lMatList,  mVecList,  rArr, sigmaRM,  sigmaERM, thetaVec, alphaRM, pMat, inverseP, decompositionH, identity, vcvMatDetArr,  negativeTwoAplusLDetArr);
        RealMatrix l0Mat = lMatList.get(tree.getRoot().getNr());
        RealVector m0Vec = mVecList.get(tree.getRoot().getNr());
        double r0 = rArr[tree.getRoot().getNr()];

        System.out.println("Display LD:");
        GeneralUtils.displayRealMatrix(lMatList.get(3));
        System.out.println("Display LE:");
        GeneralUtils.displayRealMatrix(l0Mat);
        System.out.println("Display mD:");
        GeneralUtils.displayRealVector(mVecList.get(3));
        System.out.println("Display mE:");
        GeneralUtils.displayRealVector(m0Vec);
        System.out.println("Display rD = " + rArr[3]);
        System.out.println("Display rE = " + r0);

        double loglik2 = l0Mat.preMultiply(rootValuesVec).dotProduct(rootValuesVec) + rootValuesVec.dotProduct(m0Vec) + r0;
        System.out.println("Log likelihood2 = " + loglik2);
        // expected: -7.58111239313721

        // TEST3: calculateLogP()
        pcm.initByName("tree", tree, "traits", traitValues, "alpha", alpha, "theta", theta, "sigma", sigma, "sigmae", sigmae, "root", rootValues);
        double logP = pcm.calculateLogP();
        System.out.println("Log likelihood3 = " + logP);
        // expected: -7.58111239313721
    }
}
