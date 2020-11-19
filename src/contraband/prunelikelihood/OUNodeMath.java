package contraband.prunelikelihood;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import org.apache.commons.math3.linear.*;
import outercore.parameter.KeyRealParameter;
import beast.evolution.tree.Node;
import java.util.ArrayList;
import java.util.List;

public class OUNodeMath extends CalculationNode {
    final public Input<KeyRealParameter> traitsValuesInput = new Input<>("traits","Trait values at tips.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmaValuesInput = new Input<>("sigma","Trait rate matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmaeValuesInput = new Input<>("sigmae","Population error.", Input.Validate.OPTIONAL);
    final public Input<RealParameter> alphaInput = new Input<>("alpha","An array of (nTraits * (nTraits - 1) / 2) elements, representing selection strength, off-diagonal elements in Alpha matrix.", Input.Validate.REQUIRED);
    final public Input<RealParameter> thetaInput = new Input<>("theta","An array of nTraits elements, representing optimum trait values, elements in Theta vector.", Input.Validate.REQUIRED);
    final public Input<RealParameter> rootValuesInput = new Input<>("root","Trait values at the root.");

    private RealMatrix sigmaRM;
    private RealMatrix sigmaERM;
    private RealMatrix alphaRM;
    private RealVector thetaVec;
    private RealVector rootValuesVec;
    private RealMatrix identity;
    private int nTraits;
    private int nSpecies;
    private int nodeNr;

    private List<RealMatrix> lMatList;
    private List<RealVector>  mVecList;
    private double [] rArr;
    private double [] vcvMatDetArr;
    private double [] negativeTwoAplusLDetArr;

    private RealMatrix pMat;
    private RealMatrix inverseP;
    private EigenDecomposition decompositionH;

    private RealMatrix iniMatrix;
    private RealVector iniVec;

    private RealMatrix invVCVMat;
    private double varianceRMDet;
    private RealMatrix aPlusLInv;
    private double logDetVNode;

    @Override
    public void initAndValidate() {
        // collect trait information
        KeyRealParameter traitsValues = traitsValuesInput.get();
        nTraits = traitsValues.getMinorDimension1();
        nSpecies = traitsValues.getMinorDimension2();

        nodeNr = 2 * nSpecies - 1;

        sigmaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        OUPruneUtils.populateSigmaMatrix(sigmaRM, sigmaValuesInput.get().getDoubleValues());

        alphaRM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
        OUPruneUtils.populateAlphaMatrix(alphaRM, alphaInput.get().getDoubleValues());

        thetaVec = new ArrayRealVector(new double[nTraits]);
        OUPruneUtils.populateRealVector(thetaVec, thetaInput.get().getDoubleValues());

        if (sigmaeValuesInput.get() != null) {
            sigmaERM = new Array2DRowRealMatrix(new double[nTraits][nTraits]);
            OUPruneUtils.populateSigmaMatrix(sigmaERM, sigmaeValuesInput.get().getDoubleValues());
        }

        lMatList = new ArrayList<>(nodeNr);
        mVecList = new ArrayList<>(nodeNr);
        rArr = new double[nodeNr];
        rootValuesVec = new ArrayRealVector(rootValuesInput.get().getDoubleValues());
        vcvMatDetArr = new double[nodeNr];
        negativeTwoAplusLDetArr = new double[nodeNr];

        // initiate
        iniMatrix = new Array2DRowRealMatrix(new double [nTraits][nTraits]);
        iniVec = new ArrayRealVector(new double [nTraits]);
        for (int i= 0; i < nodeNr; i ++) {
            lMatList.add(i, iniMatrix);
            mVecList.add(i, iniVec);
            rArr[i] = 0;

        }

        identity = MatrixUtils.createRealIdentityMatrix(nTraits);


    }

    // getters
    public RealMatrix getIniMatrix () {return iniMatrix;}

    public RealVector getIniVec () {return iniVec;}

    public RealMatrix getInverseVarianceMatrix () {return invVCVMat;}

    public double getDetVarianceMatrix () {return varianceRMDet ;}

    public RealMatrix getLMatForNode (int nodeIdx) {return lMatList.get(nodeIdx); }

    public RealVector getMVecForNode(int nodeIdx) {return mVecList.get(nodeIdx); }

    public double getRForNode (int nodeIdx) { return rArr[nodeIdx]; }

    public RealVector getRootValuesVec () { return rootValuesVec; }

    public RealMatrix getAlphaMatrix () {return alphaRM; }

    public RealVector getThetaVec () {return thetaVec; }

    public RealMatrix getIdentityMatrix ( ) {return identity; }

    public RealMatrix getAPlusLInv() { return aPlusLInv; }

    public double getNegativeTwoAPlusLDet() { return logDetVNode; }

    // setters
    public void setLMatForNode (int nodeIdx, RealMatrix value) { lMatList.set(nodeIdx,value); }

    public void setMVecForNode (int nodeIdx, RealVector value) { mVecList.set(nodeIdx, value); }

    public void setRForNode (int nodeIdx, double value) { rArr[nodeIdx] = value; }


    public void populateAlphaMatrix() {
        OUPruneUtils.populateAlphaMatrix(alphaRM, alphaInput.get().getDoubleValues());
    }

    public void performAlphaDecompostion () {
        decompositionH = new EigenDecomposition(alphaRM);
        pMat = decompositionH.getV();
        inverseP = new LUDecomposition(pMat).getSolver().getInverse();
    }

    public void populateVarianceCovarianceMatrix(Node node){
        invVCVMat = OUPruneUtils.getInverseVarianceRMForOU(node, vcvMatDetArr,
                sigmaRM, sigmaERM,
                pMat, inverseP, decompositionH, nTraits);
        varianceRMDet = vcvMatDetArr[node.getNr()];
    }

    public void performAPlusOperations (Node aNode, RealMatrix aMatrix) {
        // (aMat + lMat).inverse
        aPlusLInv = OUPruneUtils.getInvAPlusLRM(aNode, negativeTwoAplusLDetArr, aMatrix, lMatList.get(aNode.getNr()));

        // determinant of -2 * (aMat + lMat)
       logDetVNode = negativeTwoAplusLDetArr[aNode.getNr()];
    }

}
