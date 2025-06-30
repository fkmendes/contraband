package contraband.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.*;
import beastfx.app.util.FXUtils;
import contraband.math.NodeMath;
import javafx.geometry.Insets;
import javafx.scene.control.CheckBox;
import javafx.scene.control.TextField;
import javafx.scene.layout.HBox;


public class BMModelInputEditor extends BEASTObjectInputEditor {



	ParameterInputEditor traitsEditor;
	TextField traitsEntry;


	// vars for dealing with mean-rate delta exchange operator
	CheckBox fixMeanRatesCheckBox;

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

		pane = FXUtils.newVBox();
		pane.setPadding(new Insets(0, 5, 5, 0));
		HBox itemBox = FXUtils.newHBox();
		itemBox.setPadding(new Insets(0, 5, 5, 5));

		super.init(input, beastObject, itemNr, isExpandOption, addButtons);

	}


	//public InputEditor createTraitsEditor() {
		//NodeMath nodeMath = ((NodeMath) m_input.get());

        //final Input<?> input = nodeMath.traitsValuesInput;
		//traitsEditor = new ParameterInputEditor(doc);

		//traitsEditor.init(input, nodeMath, -1, ExpandOption.TRUE, true);
		//traitsEntry = traitsEditor.getEntry();
		//traitsEntry.setOnKeyReleased(e -> processEntry2());

		//traitsEditor.validateInput();
		//return traitsEditor;
	//}

	//void processEntry2() {
		//String[] traitValues = new String[]{traitsEntry.getText()};
		//try {
			//Double categoryCount = Double.parseDouble(traitValues[0]);
			//NodeMath s = (NodeMath) m_input.get();
			//s.getInput("traits").setValue(categoryCount, s);
			//repaint();
		//} catch (java.lang.NumberFormatException e) {
			// ignore.
		//}
	//}
}