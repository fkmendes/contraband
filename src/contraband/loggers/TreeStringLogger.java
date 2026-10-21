package contraband.loggers;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Tree;

import java.io.PrintStream;

@Description("Logs the tree as a Newick string column inside a plain tab-separated log file.")
public class TreeStringLogger extends BEASTObject implements Loggable {

    final public Input<Tree> treeInput = new Input<>("tree",
            "Tree to log as a Newick string.", Input.Validate.REQUIRED);

    private Tree tree;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
    }

    @Override
    public void init(PrintStream out) {
        out.print("tree\t");
    }

    @Override
    public void log(long sampleNr, PrintStream out) {
        out.print(tree.getRoot().toNewick() + "\t");
    }

    @Override
    public void close(PrintStream out) {
        // nothing to close
    }
}
