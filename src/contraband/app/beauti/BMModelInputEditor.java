package contraband.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.*;
import beastfx.app.util.FXUtils;
import contraband.math.NodeMath;
import javafx.geometry.Insets;
import javafx.scene.control.CheckBox;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;



public class BMModelInputEditor extends BEASTObjectInputEditor {


	ParameterInputEditor sigmasqEditor;
	TextField sigmasqEntry;
	protected SmallLabel sigmasqLabel;

	// vars for dealing with mean-rate delta exchange operator
	CheckBox ratesCheckBox;

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


	public InputEditor createSigmasqEditor() {
		NodeMath nodeMath = ((NodeMath) m_input.get());

		final Input<?> input = nodeMath.sigmasqInput;
		sigmasqEditor = new ParameterInputEditor(doc);
		sigmasqEditor.init(input, nodeMath, -1, ExpandOption.FALSE, true);
		sigmasqEntry = sigmasqEditor.getEntry();

		sigmasqEditor.validateInput();
		return sigmasqEditor;
	}


}