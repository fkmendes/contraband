package test;

import beast.util.TreeParser;

public class OUOneTraitAnalyticalTestDriver {

	public static void main(String[] args) {
		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=1]:2.0)[&Regime=1]:1.0, sp4[&Regime=1]:3.0)[&Regime=1];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
	}

}
