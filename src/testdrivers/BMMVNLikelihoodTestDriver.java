package testdrivers;

import java.util.Arrays;
import java.util.List;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import contraband.BMMVNLikelihoodOneTrait;
import contraband.OneValueContTraits;

public class BMMVNLikelihoodTestDriver {

	public static void main(String[] args) {
		// tree
		String treeStr = "(((sp1[&Regime=1]:1.0, sp2[&Regime=1]:1.0)[&Regime=1]:1.0, sp3[&Regime=2]:2.0)[&Regime=0]:1.0, sp4[&Regime=0]:3.0)[&Regime=0];";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		// initializing data
		String[] spNames = new String[] { "sp1", "sp2", "sp3", "sp4" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		String oneTraitValues = "sp1=4.1,sp2= 4.5,sp3=5.9,sp4 =0.0 ";
		OneValueContTraits oneTraitData = new OneValueContTraits();
		oneTraitData.initByName("nTraits", 1, "taxa", taxonSet, "traitValues", oneTraitValues);
		
		// sigmasq
		Double[] sigmasqInput = new Double[] { 1.4822794118 };
		RealParameter sigmasq = new RealParameter(sigmasqInput);
		
		// mean vector
		Double[] meanVectorInput = new Double[] { 3.079142, 3.079142, 3.079142, 3.079142 };
		RealParameter mean = new RealParameter(meanVectorInput);
		
		// likelihood
		BMMVNLikelihoodOneTrait BMLk = new BMMVNLikelihoodOneTrait();
		BMLk.initByName("tree", myTree, "sigmasq", sigmasq, "mean", mean, "data", oneTraitData);
		System.out.println(BMLk.calculateLogP());
		}
}
