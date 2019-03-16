package test;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Before;
import org.junit.Test;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import contraband.OneValueContTraits;

public class BMMVNLikelihoodOneTraitTest {

	@Before
	public void setUp() throws Exception {
		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		List<Node> allLeafNodes = myTree.getRoot().getAllLeafNodes();
		
		// initializing data
		String[] spNames = new String[] { "sp1", "sp2", "sp3", "sp4" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		String oneTraitValues = "sp1=4.1,sp2= 4.5,sp3=5.9,sp4 =0.0 ";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "taxa", taxonSet, "traitValues", oneTraitValues);
		
		
		
		
	}

	@Test
	public void test() {
		fail("Not yet implemented");
	}

}
