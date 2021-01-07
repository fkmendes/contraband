package testdrivers;

import java.util.Arrays;

public class LiabilityLikelihoodTestDriver {

    public static void main(String[] args) {
        double[] liabilities1 = new double[]{
                4.0, 5.0,
                9.0, 10.0,
                14.0, 15.0};
        System.out.println("liability1 = " + "\n" + Arrays.toString(liabilities1));

        double[] liabilities2 = new double[]{
                16.0, 17.0,
                18.0, 19.0,
                20.0, 21.0};
        System.out.println("liability2 = " + "\n" + Arrays.toString(liabilities2));

        double[] traitValuesArr = new double[]{
                1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0,
                6.0, 7.0, 8.0, 0.0, 0.0, 0.0, 0.0,
                11.0, 12.0, 13.0, 0.0, 0.0, 0.0, 0.0};
        System.out.println("initial = " + "\n" +  Arrays.toString(traitValuesArr));

        int traitNr = 2;
        int totalNr = 7;

        int index = 3;
        for (int i = 0; i < 3; i++) {
            System.arraycopy(liabilities1, i * traitNr, traitValuesArr, i * totalNr + index, traitNr);
        }
        System.out.println("populated1 = " + "\n" + Arrays.toString(traitValuesArr));

        index = 5;
        for (int i = 0; i < 3; i++) {
            System.arraycopy(liabilities2, i * traitNr, traitValuesArr, i * totalNr + index, traitNr);
        }
        System.out.println("populated2 = " + "\n" + Arrays.toString(traitValuesArr));
    }

}
