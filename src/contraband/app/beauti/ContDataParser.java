package contraband.app.beauti;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.datatype.StandardData;
import beast.base.inference.parameter.RealParameter;
import beast.base.parser.NexusParser;
import beast.base.parser.XMLProducer;

import java.io.*;
import java.util.*;

public class ContDataParser {


    protected Integer taxaNr;
    protected Integer charNr;

    public Alignment m_alignment;
    public List<String> taxa;
    List<Taxon> taxonList = new ArrayList<>();
    static Set<String> g_sequenceIDs;

    static {
        g_sequenceIDs = new HashSet<>();
    }

    public void parseFile(final File file) throws IOException {
        final String fileName = file.getName().replaceAll(".*[\\/\\\\]", "").replaceAll("\\..*", "");

        parseFile(fileName, new FileReader(file));
    }

    public void parseFile(final String id, final Reader reader) throws IOException {

        final BufferedReader fin;
        if (reader instanceof BufferedReader) {
            fin = (BufferedReader) reader;
        } else {
            fin = new BufferedReader(reader);
        }
        String firstLine = fin.readLine();
        String[] firstLineValues = firstLine.split("\t");
        taxaNr = Integer.parseInt(firstLineValues[0].split("=")[1]);
        charNr = Integer.parseInt(firstLineValues[1].split("=")[1]);


        m_alignment = parseDataBlock(fin);
        m_alignment.setID(id);
    } // parseFile



    private boolean taxonListContains(String taxon) {
        for (Taxon t : taxonList) {
            if (t.getID().equals(taxon)) {
                return true;
            }
        }
        return false;
    }

    /**
     * parse data block and create Alignment *
     */
    public Alignment parseDataBlock(final BufferedReader fin) throws IOException {
        taxa = new ArrayList<>();
        final Alignment alignment = new Alignment();
        alignment.dataTypeInput.setValue("user defined", alignment);
        ContinuousData type = new ContinuousData();
        type.setInputValue("nrOfStates", charNr);
        type.initAndValidate();
        alignment.setInputValue("userDataType", type);
        String ambiguitiesStr = "";

        for (int i = 0; i < taxaNr; i++) {
            String line = fin.readLine();



            String[] strValues = line.split("\t");
            String taxon = strValues[0];
            if (!taxa.contains(taxon)) {
                taxa.add(taxon);
            }

            if (!taxonListContains(taxon)) {
                taxonList.add(new Taxon(taxon));
            }

            int idx = 0;
            final StringBuilder traitValues = new StringBuilder();
            for (int j = 1; j < charNr + 1; j++) {
                if (idx == 0) {
                    traitValues.append("{").append(strValues[j]).append("}");
                } else {
                    traitValues.append(" {").append(strValues[j]).append("}");
                }
                ambiguitiesStr += strValues[j]+ " ";
                idx = idx + 1;
            }

            String data = traitValues.toString();
            final Sequence sequence = new Sequence();
            sequence.init(charNr, taxon, data);
            sequence.setID(generateSequenceID(taxon));
            alignment.sequenceInput.setValue(sequence, alignment);
        }

        alignment.userDataTypeInput.get().initByName("ambiguities", ambiguitiesStr);
        alignment.initAndValidate();
        return alignment;
    }



    public static void main(final String[] args) {
        try {
            final NexusParser parser = new NexusParser();
            parser.parseFile(new File(args[0]));
            if (parser.taxa != null) {
                System.out.println(parser.taxa.size() + " taxa");
                System.out.println(Arrays.toString(parser.taxa.toArray(new String[parser.taxa.size()])));
            }
            if (parser.trees != null) {
                System.out.println(parser.trees.size() + " trees");
            }
            if (parser.m_alignment != null) {
                final String xml = new XMLProducer().toXML(parser.m_alignment);
                System.out.println(xml);
            }
            if (parser.traitSet != null) {
                final String xml = new XMLProducer().toXML(parser.traitSet);
                System.out.println(xml);
            }
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    } // main

    public static String generateSequenceID(final String taxon) {
        String id = "seq_" + taxon;
        int i = 0;
        while (g_sequenceIDs.contains(id + (i > 0 ? i : ""))) {
            i++;
        }
        id = id + (i > 0 ? i : "");
        g_sequenceIDs.add(id);
        return id;
    }
}
