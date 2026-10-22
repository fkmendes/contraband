package contraband.loggers;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import contraband.math.NodeMath;

import java.io.PrintStream;

@Description("Logs body mass (trait mean given all the data, from NodeMath) and branch rate " +
        "for every node into a tab-separated log file. " +
        "Use TreeStringLogger in the same <logger> block to also record the tree.")
public class InternalNodesBodyMassLogger extends BEASTObject implements Loggable {

    final public Input<NodeMath> nodeMathInput = new Input<>("nodeMath",
            "NodeMath holding internal-node trait expectations.", Input.Validate.REQUIRED);

    final public Input<BranchRateModel.Base> clockModelInput = new Input<>("branchRateModel",
            "Branch rate model used to compute branch rates.", Input.Validate.REQUIRED);

    final public Input<Tree> treeInput = new Input<>("tree",
            "Tree whose nodes will be logged.", Input.Validate.REQUIRED);

    final public Input<Integer> traitIndexInput = new Input<>("traitIndex",
            "0-based index of the trait to log as body mass. Default: 0.", 0);

    private NodeMath nodeMath;
    private BranchRateModel.Base clockModel;
    private Tree tree;
    private int nSpecies;
    private int traitIndex;
    // true if any node in the tree has a non-empty ID string (taxa name or
    // internal node label set from a labelled Newick via IsLabelledNewick="true").
    private boolean hasLabels;

    @Override
    public void initAndValidate() {
        nodeMath   = nodeMathInput.get();
        clockModel = clockModelInput.get();
        tree       = treeInput.get();
        traitIndex = traitIndexInput.get();
        nSpecies   = tree.getLeafNodeCount();

        // Scan once: if any node carries a non-empty ID, use getID() for all;
        // otherwise fall back to getNr() for all.
        hasLabels = false;
        for (Node node : tree.getNodesAsArray()) {
            String id = node.getID();
            if (id != null && !id.isEmpty()) {
                hasLabels = true;
                break;
            }
        }
    }

    /**
     * Returns a label for this node.
     * If the tree has no labels at all: "node_<nr>" for every node.
     * If the tree has labels (taxa names and/or internal node labels):
     *   use getID(); fall back to "node_<nr>" for any node whose ID is unset.
     */
    private String nodeLabel(Node node) {
        if (!hasLabels) {
            return "node_" + node.getNr();
        }
        String id = node.getID();
        return (id != null && !id.isEmpty()) ? id : "node_" + node.getNr();
    }

    @Override
    public void init(PrintStream out) {
        int nodeCount = tree.getNodeCount();

        // column headers — one bodyMass and one rate column per node, no tree column
        for (int i = 0; i < nodeCount; i++) {
            String label = nodeLabel(tree.getNode(i));
            out.print("bodyMass_" + label + "\t");
            out.print("rate_"     + label + "\t");
        }

    }

    @Override
    public void log(long sampleNr, PrintStream out) {
        int nodeCount = tree.getNodeCount();
        nodeMath = nodeMathInput.get();

        for (int i = 0; i < nodeCount; i++) {
            Node node = tree.getNode(i);

            // mean given all the data (tipward pass); for tips with data, the observed value
            double bodyMass = nodeMath.getConditionalMean(node.getNr(), traitIndex);
            double rate = clockModel.getRateForBranch(node);

            out.print(bodyMass + "\t");
            out.print(rate     + "\t");
        }

    }

    @Override
    public void close(PrintStream out) {
        // nothing to close
    }
}
