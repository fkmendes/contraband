package contraband.app.beauti;


import java.io.File;
import java.util.*;


import beast.base.core.Description;
import beast.base.core.ProgramStatus;
import beastfx.app.inputeditor.BeautiAlignmentProvider;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.Alert;
import beastfx.app.util.ExtensionFileFilter;
import beastfx.app.util.FXUtils;
import javafx.scene.control.ButtonType;
import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.FilteredAlignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.StandardData;
import beast.base.evolution.datatype.UserDataType;
import beast.base.parser.NexusParser;
import beast.base.parser.PartitionContext;


@Description("Class for creating new partitions for morphological data to be edited by AlignmentListInputEditor")
public class BeautiContrabandAlignmentProvider extends BeautiAlignmentProvider {

	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc) {
		String[] exts = { "nex", "nxs", "nexus"};
		File [] files = FXUtils.getLoadFiles("Load Alignment File",
				new File(ProgramStatus.g_sDir), "Alignment files", exts);

		if (files != null && files.length > 0) {

			List<BEASTInterface> oneTraitData = getAlignments(doc, files);
			return oneTraitData;

		}
		return null;
	}

	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {
		if (files == null) {
			// merge "+ button" and "drag drop" function
			return getAlignments(doc);
		}
		List<BEASTInterface> selectedPlugins = new ArrayList<BEASTInterface>();
		for (File file : files) {
			String fileName = file.getName();

				ContDataParser parser = new ContDataParser();

				try {
					parser.parseFile(file);
					selectedPlugins.add(parser.m_alignment);
				} catch (Exception ex) {
					ex.printStackTrace();
					Alert.showMessageDialog(null, "Loading of " + fileName + " failed: " + ex.getMessage());
					return null;
				}

		}


		return selectedPlugins;

	}



	public void processAlignment(Alignment alignment, List<BEASTInterface> filteredAlignments, boolean ascertained, BeautiDoc doc) throws Exception {

		int nrOfStates = 3;
		String tree = alignment.getID();
		String clock = alignment.getID();
		String ID = alignment.getID() + nrOfStates;


			// create data type
			DataType.Base dataType = null;
			if (alignment.getDataType() instanceof StandardData) {
				// determine state space size by interrogating StandardData
				// data-type
				StandardData base = (StandardData) alignment.getDataType();
				dataType = new StandardData();
				((StandardData) dataType).initByName("nrOfStates", nrOfStates,
						"ambiguities", base.listOfAmbiguitiesInput.get());
			}

			String name = alignment.getID() + nrOfStates;
			dataType.setID("contDataType." + name);
			doc.addPlugin(dataType);

		FilteredAlignment data = new FilteredAlignment();
		data.initByName("data", alignment, "filter", "1", "userDataType", dataType);
		data.setID(ID);
		doc.addPlugin(data);


			// link trees and clock models
			PartitionContext context = new PartitionContext(name, name, clock, tree);

			// create treelikelihood for each state space
			try {
				doc.addAlignmentWithSubnet(context, template.get());
//				GeneralSubstitutionModel smodel = (GeneralSubstitutionModel) doc.pluginmap.get("morphSubstModel.s:" + name);
//				((RealParameter) smodel.ratesInput.get()).setDimension(nrOfStates * (nrOfStates - 1)/2);
//				smodel.frequenciesInput.get().frequenciesInput.get().setDimension(nrOfStates);
//				SiteModelInterface.Base sitemodel = (SiteModelInterface.Base) doc.pluginmap.get("morphSiteModel.s:" + name);
//				sitemodel.substModelInput.setValue(smodel, sitemodel);
			} catch (Exception e) {
				e.printStackTrace();
			}

			filteredAlignments.add(alignment);


	}


	@Override
	public int matches(Alignment alignment) {
		if (alignment.userDataTypeInput.get() != null && alignment.userDataTypeInput.get() instanceof StandardData) {
			return 20;
		}
		return 0;
	}


}
