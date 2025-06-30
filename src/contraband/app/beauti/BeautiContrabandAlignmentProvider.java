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

		List<BEASTInterface> processedAlignments = new ArrayList<>();
		for (BEASTInterface o : selectedPlugins) {
			try {
			if (o instanceof Alignment) {

				processAlignment((Alignment) o, processedAlignments, doc);
			}
			} catch (Exception e) {
				Alert.showMessageDialog(null, "Something went wrong converting the alignment: " + e.getMessage());
				e.printStackTrace();
				return null;
			}
		}



		return selectedPlugins;

	}



	public void processAlignment(Alignment alignment, List<BEASTInterface> processedAlignments, BeautiDoc doc) throws Exception {
		Map<Integer, List<Integer>> stateSpaceMap = new HashMap<Integer, List<Integer>>();


		ContinuousData dataType = (ContinuousData) alignment.getDataType();

		String name = alignment.getID() ;
		String tree = alignment.getID() ;
		String clock = alignment.getID() ;

		doc.addPlugin(dataType);
		doc.addPlugin(alignment);
		// link trees and clock models
		PartitionContext context = new PartitionContext(name, name, clock, tree);

		try {
			doc.addAlignmentWithSubnet(context, template.get());
		} catch (Exception e) {
			e.printStackTrace();
		}

		processedAlignments.add(alignment);

	}


	@Override
	public int matches(Alignment alignment) {
		if (alignment.userDataTypeInput.get() != null && alignment.userDataTypeInput.get() instanceof ContinuousData) {
			return 20;
		}
		return 0;
	}


}
