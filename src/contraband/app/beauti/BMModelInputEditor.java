package contraband.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.*;
import contraband.math.NodeMath;
import contraband.prunelikelihood.BMPruneLikelihood;
import contraband.prunelikelihood.BMPruneShrinkageLikelihood;
import javafx.scene.control.TextField;


public class BMModelInputEditor extends BEASTObjectInputEditor {


	BEASTObjectInputEditor nodeMathEditor;
	ParameterInputEditor popTraitsEditor;
	TextField popTraitsEntry;

	public BMModelInputEditor() {
		super();
	}
	public BMModelInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return NodeMath.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
					 ExpandOption isExpandOption, boolean addButtons) {
		m_bAddButtons = addButtons;
		m_input = input;
		m_beastObject = beastObject;

		super.init(input, beastObject, itemNr, isExpandOption, addButtons);
	}


	public InputEditor createPopulationTraitsEditor() {

		NodeMath nodeMath = (NodeMath) m_input.get();

		final Input<?> input = nodeMath.populationTraitsInput;

		popTraitsEditor = new ParameterInputEditor(doc);
		popTraitsEditor.init(input, nodeMath, -1, ExpandOption.TRUE, true);
		popTraitsEntry = popTraitsEditor.getEntry();

		popTraitsEditor.validateInput();
		return popTraitsEditor;
	}


}