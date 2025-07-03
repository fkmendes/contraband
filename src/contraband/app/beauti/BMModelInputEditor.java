package contraband.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.*;
import contraband.prunelikelihood.BMPruneLikelihood;


public class BMModelInputEditor extends BEASTObjectInputEditor {


	BEASTObjectInputEditor nodeMathEditor;


	public BMModelInputEditor() {
		super();
	}
	public BMModelInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return BMPruneLikelihood.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
					 ExpandOption isExpandOption, boolean addButtons) {
		m_bAddButtons = addButtons;
		m_input = input;
		m_beastObject = beastObject;

		super.init(input, beastObject, itemNr, isExpandOption, addButtons);
	}


	public InputEditor createNodeMathEditor() {
		BMPruneLikelihood lik = (BMPruneLikelihood) m_input.get();
		final Input<?> input = lik.nodeMathInput;

		nodeMathEditor = new BEASTObjectInputEditor(doc);
		nodeMathEditor.init(input, lik, -1, ExpandOption.TRUE, true);

		nodeMathEditor.validateInput();
		return nodeMathEditor;
	}


}