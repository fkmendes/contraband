package testdrivers;

import beast.util.TreeParser;
import contraband.GeneralUtils;
import contraband.PruneUtils;
import org.apache.commons.math3.linear.*;
import beast.evolution.tree.Node;
import java.util.ArrayList;
import java.util.List;

public class PruneUtilsTestDriver {
    public static void main(String[] args) {
        // initialize bifurcating tree
        // set branch lengths to 1.0 for simplicity
        // tree
        String treeStr =  "((t3:1.0,t4:1.0):1.0,(t1:1.0,t2:1.0):1.0):0.0;";
        TreeParser tree = new TreeParser(treeStr, false, false, true, 0);

        double[][] evolRateMat2DArray = {{9.404673, -9.088870, 7.848381,  9.1366669},
                                         {0.000000,  5.281055, 1.028700, -0.9333169},
                                         {0.000000,  0.000000, 4.566147,  3.5514127},
                                         {0.000000,  0.000000, 0.000000,  5.7263340}}; // use a simple example for four traits
        RealMatrix evolRateMat = new Array2DRowRealMatrix(evolRateMat2DArray);

        // initiate lists
        RealMatrix initRM = new Array2DRowRealMatrix(new double[4][4]);
        RealVector iniVec = new ArrayRealVector(new double[4]);
        List<RealMatrix> vcvMatList = new ArrayList<>(tree.getNodeCount());
        List<RealMatrix> aMatList = new ArrayList<>(tree.getNodeCount());
        List<RealMatrix> eMatList = new ArrayList<>(tree.getNodeCount());
        List<RealMatrix> cMatList = new ArrayList<>(tree.getNodeCount());
        List<RealVector> dVecList = new ArrayList<>(tree.getNodeCount());
        double[] fArrary = new double[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            vcvMatList.add(i, new Array2DRowRealMatrix(new double[10][10]));
            aMatList.add(i, initRM);
            eMatList.add(i, initRM);
            cMatList.add(i, initRM);
            dVecList.add(i, iniVec);
        }

        // the node to calculate
        Node aNode = tree.getNode(1);

        // block for matrix VCV
        System.out.println("Printing matrix VCVMat at node 1:");
        PruneUtils.setVCVMatForBranchBM(aNode, evolRateMat, vcvMatList);
        RealMatrix vCVMat = vcvMatList.get(aNode.getNr());
        GeneralUtils.displayRealMatrix(vCVMat);

        RealMatrix inverseVCVMat = new Array2DRowRealMatrix(new double [4][4]);
        try {
            LUDecomposition VMatLUD = new LUDecomposition(vCVMat);
            inverseVCVMat = VMatLUD.getSolver().getInverse();
        }
        catch (SingularMatrixException e) {
            System.out.println("vCVMat at this node is singular");
        }

        // block for matrix A
        System.out.println("Printing matrix A:");
        // write code to initialize aMat from PruneUtils
        PruneUtils.setAMat(aNode, vCVMat, aMatList);
        RealMatrix aMat = aMatList.get(aNode.getNr());
        GeneralUtils.displayRealMatrix(aMat);

        // block for matrix E
        System.out.println("Printing matrix E:");
        // write code to initialize eMat from PruneUtils
        RealMatrix phiMat = MatrixUtils.createRealIdentityMatrix(4);
        PruneUtils.setEMatOU(aNode, inverseVCVMat, phiMat.transpose(), eMatList);
        RealMatrix eMat = eMatList.get(aNode.getNr());
        GeneralUtils.displayRealMatrix(eMat);

        // block for matrix C
        System.out.println("Printing matrix C:");
        // write code to initialize cMat from PruneUtils
        PruneUtils.setCMatOU(aNode, eMat, phiMat, cMatList);
        RealMatrix cMat = cMatList.get(aNode.getNr());
        GeneralUtils.displayRealMatrix(cMat);

        // block for vector B
        System.out.println("Printing vector B:");
        // write code to initialize bVec
        RealVector bVec = new ArrayRealVector(new double [4]);
        GeneralUtils.displayRealVector(bVec);

        // block for vector D
        System.out.println("Printing vector D:");
        // write code to initialize dVec from PruneUtils
        RealVector omegaVec = new ArrayRealVector(new double [4]);
        PruneUtils.setDVecOU(aNode, eMat, omegaVec, dVecList);
        RealVector dVec = dVecList.get(aNode.getNr());
        GeneralUtils.displayRealVector(dVec);

        // block for double f
        double vCVMatDet = 0.0;
        try {
            LUDecomposition VMatLUD = new LUDecomposition(vCVMat);
            vCVMatDet = VMatLUD.getDeterminant();
        }
        catch (SingularMatrixException e) {
            System.out.println("vCVMat at this node is singular");
        }
        PruneUtils.setF(aNode,4, vCVMatDet, fArrary);
        double f = fArrary[aNode.getNr()];
        System.out.println("Printing double f:" + f);

    }
}
