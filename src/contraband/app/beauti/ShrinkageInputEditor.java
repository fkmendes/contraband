package contraband.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.*;
import beastfx.app.util.FXUtils;
import contraband.math.NodeMath;
import contraband.prunelikelihood.BMPruneShrinkageLikelihood;
import javafx.geometry.Insets;
import javafx.scene.control.TextField;
import javafx.scene.layout.HBox;


public class ShrinkageInputEditor extends BEASTObjectInputEditor {


	ParameterInputEditor popTraitsEditor;
	TextField popTraitsEntry;

	ParameterInputEditor sigmasqEditor;
	TextField sigmasqEntry;

	DoubleInputEditor popVarEditor;
	TextField popVarEntry;

	public ShrinkageInputEditor() {
		super();
	}
	public ShrinkageInputEditor(BeautiDoc doc) {
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


	public InputEditor createPopulationTraitsEditor() {

		BMPruneShrinkageLikelihood lik = (BMPruneShrinkageLikelihood) m_input.get();

		final Input<?> input = lik.populationTraitsInput;

		popTraitsEditor = new ParameterInputEditor(doc);
		popTraitsEditor.init(input, lik, -1, ExpandOption.TRUE, true);
		popTraitsEntry = popTraitsEditor.getEntry();

		popTraitsEditor.validateInput();
		return popTraitsEditor;
	}

	/*
	public InputEditor createPopVarEditor() {
		BMPruneShrinkageLikelihood lik = (BMPruneShrinkageLikelihood) m_input.get();
		final Input<?> input = lik.popVarInput;
		popVarEditor = new DoubleInputEditor(doc);

		popVarEditor.init(input, lik, -1, ExpandOption.FALSE, true);
		popVarEntry = popVarEditor.getEntry();
		return popVarEditor;
	}
*/

	public InputEditor createNodeMathEditor() {
		BMPruneShrinkageLikelihood lik = (BMPruneShrinkageLikelihood) m_input.get();
		final Input<?> input = lik.nodeMathInput;

		NodeMath nodeMath = ((NodeMath) input.get());

		final Input<?> sigmasqInput = nodeMath.sigmasqInput;
		sigmasqEditor = new ParameterInputEditor(doc);
		sigmasqEditor.init(sigmasqInput, nodeMath, -1, ExpandOption.FALSE, true);
		sigmasqEntry = sigmasqEditor.getEntry();

		sigmasqEditor.validateInput();
		return sigmasqEditor;

	}

}