package contraband.app.beauti;

import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ContinuousData extends DataType.Base {
    final public Input<Integer> maxNrOfStatesInput = new Input<>("nrOfStates", "specifies the maximum number of " +
            "character states in data matrix or in the filtered alignment");
    final public Input<String> listOfAmbiguitiesInput = new Input<>("ambiguities", "all possible ambiguities presented " +
            "as space separated sets of ordered elements. Elements are digits 0..9.");

    private String[] ambiguities = {};
    private ArrayList<String> codeMapping;

    private int ambCount;

    @Override
    public void initAndValidate() {
        stateCount = maxNrOfStatesInput.get();

        mapCodeToStateSet = null;
        codeLength = -1;
        codeMap = null;
        createCodeMapping();
    }

    private void createCodeMapping() {
        if (listOfAmbiguitiesInput.get() != null) {
            ambiguities = listOfAmbiguitiesInput.get().split(" ");


            ambCount = ambiguities.length;
            codeMapping = new ArrayList<>();
            for (int i=0; i< ambCount; i++) {
                codeMapping.add(ambiguities[i]);
            }

            mapCodeToStateSet = new int[ambCount][];
            for (int i = 0; i < ambCount; i++) {
                int[] stateSet = new int[codeMapping.get(i).length()];
                for (int k = 0; k < stateSet.length; k++) {
                    stateSet[k] = (codeMapping.get(i).charAt(k) - '0');
                }
                mapCodeToStateSet[i] = stateSet;
            }

            // TODO: is this the correct way to deal with stateCount == -1?
            int n = stateCount >= 0 ? stateCount : 10;
            int[] stateSet = new int[n];
            for (int i = 0; i < n; i++) {
                stateSet[i] = i;
            }
            // GAP_CHAR
            mapCodeToStateSet[mapCodeToStateSet.length - 2] = stateSet;
            // MISSING_CHAR
            mapCodeToStateSet[mapCodeToStateSet.length - 1] = stateSet;
        }
    }

    @Override
    public int[] getStatesForCode(int state) {
        if (state >= 0) {
            return mapCodeToStateSet[state];
        } else {
            return mapCodeToStateSet[mapCodeToStateSet.length - 1];
        }
    }


    @Override
    public List<Integer> stringToEncoding(String data) {
        List<Integer> sequence;
        sequence = new ArrayList<>();
        // remove spaces
        data = data.replaceAll("\\s", "");

        ArrayList<Integer> amb = new ArrayList<>();
        boolean readingAmb=false;
        for (byte c : data.getBytes()) {
            if (!readingAmb) {
                switch (c) {
                    case GAP_CHAR:
                    case MISSING_CHAR:
                        String missing = Character.toString(MISSING_CHAR);
                        sequence.add(codeMapping.indexOf(missing));
                        break;
                    case '{':
                        readingAmb = true;
                        amb.clear();
                        break;
                    default:
                        sequence.add(Integer.parseInt((char) c + ""));
                }
            } else {
                if (c != '}') {
                    amb.add(Integer.parseInt((char) c + "") );
                } else {
                    readingAmb = false;
                    String ambStr = "";
                    for (Integer a : amb) {
                        ambStr += Integer.toString(a);
                    }
                    int x = codeMapping.indexOf(ambStr);
                    sequence.add(codeMapping.indexOf(ambStr));
                }

            }

        }

        return sequence;
    }


    @Override
    public String getTypeDescription() {
        return "continuous";
    }

    @Override
    public String getCharacter(int state) {
        return codeMapping.get(state);
    }


    @Override
    public boolean hasConstantCodeLength() {
        return false;
    }
}
