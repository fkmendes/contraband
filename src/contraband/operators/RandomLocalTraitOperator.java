/*
* File BitFlipOperator.java
*
* Copyright (C) 2010 Joseph Heled jheled@gmail.com
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/

package contraband.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import contraband.clock.RandomLocalColorModel;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */

@Description("Selects a trait to operate on with a random walk. With probability p it will choose a trait associated " +
        "with an indicator of 1, otherwise it will choose a trait associated with an indicator of 0.")
public class RandomLocalTraitOperator extends Operator {

	final public Input<BooleanParameter> shiftIndicatorsInput = new Input<>("shiftIndicators", "the boolean indicator parameter, with one boolean per node (without root).", Validate.REQUIRED);
	
	final public Input<RealParameter> traitsInput = new Input<>("traits", "the trait (rate, optimum, etc.) associated with each tree node (including root).", Validate.REQUIRED);
			
    final public Input<RandomLocalColorModel> randomColorModelInput =
            new Input<>("randomColorModel",
                    "the random color model on whose parameters to operate.",
                    Input.Validate.REQUIRED);

    final public Input<Double> pInput = new Input<>("p", "The probability of choosing a trait with an associated 1 (including the root trait when applicable).", 0.5, Validate.OPTIONAL);

    final public Input<Double> windowSizeInput =
            new Input<>("windowSize", "the size of the window both up and down.", Input.Validate.REQUIRED);
    

    List<Integer> onPositions = new ArrayList<>();
    List<Integer> offPositions = new ArrayList<>();

    RandomLocalColorModel rlcModel;
    BooleanParameter indicatorParam;
    RealParameter traitParam;
    double windowSize = 1;

    @Override
	public void initAndValidate() {
        rlcModel = randomColorModelInput.get();
        // indicatorParam = rlcModel.indicatorParamInput.get();
        // traitParam = rlcModel.rateParamInput.get();
    	
    	indicatorParam = shiftIndicatorsInput.get();
    	traitParam = traitsInput.get();
    	
        windowSize = windowSizeInput.get();
    }

    /**
     * A bit move picks a random 'on' bit and moves it to a new position that was previously 'off'. The original position becomes 'off'.
     * This is effectively two bit flips, one from 'on' to 'off' and one from 'off' to 'on'.
     */

    @Override
    public double proposal() {
        populateOnOffPositions();

        if (Randomizer.nextDouble() < pInput.get()) {
            // operate on an 'on' trait

            int index = Randomizer.nextInt(onPositions.size());
            return proposeNewTrait(onPositions.get(index));

        } else {
            // operator on an 'off' trait
            int index = Randomizer.nextInt(offPositions.size());
            return proposeNewTrait(offPositions.get(index));
        }
    }

    private double proposeNewTrait(int nodeNumber) {

        double value = traitParam.getValue(nodeNumber);
        double newValue = value + Randomizer.nextDouble() * 2 * windowSize - windowSize;

        if (newValue < traitParam.getLower() || newValue > traitParam.getUpper()) {
            return Double.NEGATIVE_INFINITY;
        }

        traitParam.setValue(nodeNumber, newValue);

        return 0.0;

    }

    private void populateOnOffPositions() {
        onPositions.clear();
        offPositions.clear();

        if (rlcModel.includeRootInput.get()) {
            onPositions.add(rlcModel.treeInput.get().getRoot().getNr());
        }
        for (int i = 0; i < indicatorParam.getDimension(); i++) {
            if (indicatorParam.getValue(i)) {
                onPositions.add(i);
            } else {
                offPositions.add(i);
            }
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return windowSize;
    }

    @Override
    public void setCoercableParameterValue(double value) {
        windowSize = value;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */
    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
        double delta = calcDelta(logAlpha);

        delta += Math.log(windowSize);
        windowSize = Math.exp(delta);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = windowSize * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }
}

