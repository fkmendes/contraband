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

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Logger;
import beast.base.inference.Operator;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Alexei Drummond
 */

@Description("Move k 'on' bits in an array of boolean bits. The number of 'on' bits remains constant under this operation.")
public class BitMoveOperator extends Operator {
    
	final public Input<Integer> kInput = new Input<>("k", "the number of 'on' bits to shift", 1);
    final public Input<BooleanParameter> parameterInput = new Input<>("parameter", "the parameter to operate a bit move on.", Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
    final public Input<Boolean> includeSampledAncestorsInput = new Input<>("includeSampledAncestors", "true if allowed to move to sampled ancestors (i.e, nodes whose subtending branches have length zero).", false, Validate.OPTIONAL);
    
    Tree tree;
    List<Integer> onPositions = new ArrayList<>();
    List<Integer> offPositions = new ArrayList<>();

    @Override
	public void initAndValidate() {
    	tree = treeInput.get();
    }

    /**
     * A bit move picks a random 'on' bit and moves it to a new position that was previously 'off'. The original position becomes 'off'.
     * This is effectively two bit flips, one from 'on' to 'off' and one from 'off' to 'on'.
     */

    @Override
    public double proposal() {

        final BooleanParameter p = parameterInput.get();

        final int dim = p.getDimension();

        onPositions.clear();
        offPositions.clear();
        for (int i = 0; i < dim; i++) {
            if (p.getValue(i)) onPositions.add(i);
                else offPositions.add(i);
        }

        if (onPositions.size() == 0 || offPositions.size() == 0) {
            Log.warning("BitMoveOperator has no valid moves. Rejecting.");
            return Double.NEGATIVE_INFINITY;
        }

        for (int i = 0; i < kInput.get(); i++) {

        	int onPos = -1;
        	while (true) {
        		onPos = Randomizer.nextInt(onPositions.size());
        		if (!tree.getNode(onPos).isDirectAncestor()) { break; }
        	}
            Integer on = onPositions.get(onPos);

            int offPos = Randomizer.nextInt(offPositions.size());
            Integer off = offPositions.get(offPos);

            p.setValue(on,false);
            p.setValue(off,true);

            onPositions.remove(onPos);
            onPositions.add(off);
            offPositions.remove(offPos);
            offPositions.add(on);
        }

        return 0.0;
    }
}

