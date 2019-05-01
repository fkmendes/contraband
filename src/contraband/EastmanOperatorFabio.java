package contraband;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import beast.util.TreeParser;

public class EastmanOperatorFabio extends Operator {
    final public Input<TreeParser> treeInput = new Input<>("tree", "The number of colours.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> nColorInput = new Input<>("nColor", "The number of colours.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> colorAssignmentInput = new Input<>("colorAssignments", "Integers representing colors, one per branch.", Input.Validate.REQUIRED);
    
    @Override
    public void initAndValidate() {
    	;
    }

    @Override
    public double proposal() {
    	TreeParser tree = treeInput.get();
    	int nSpp = tree.getLeafNodeCount();
    	int totalNColor = tree.getNodeCount() - 1;
    	Integer[] allColors = new Integer[totalNColor];
    	for (int i=0; i<totalNColor; ++i) { 
    		allColors[i] = i;
    	}
    	Set<Integer> inactiveColorsSet = new HashSet<>(Arrays.asList(allColors));
        
    	Integer[] colorAssignments = colorAssignmentInput.get().getValues();
        Set<Integer> activeColorsSet = new HashSet<>(Arrays.asList(colorAssignments));
        inactiveColorsSet.removeAll(activeColorsSet); // hashset of inactive colors populated
        
        System.out.println("\ncolorAssignments=" + Arrays.toString(colorAssignments));
        for (int activeColor: activeColorsSet) {
        	System.out.println("Active color: " + activeColor);
        }
        for (int inactiveColor: inactiveColorsSet) {
        	System.out.println("Inactive color: " + inactiveColor);
        }
        
        // getting node to compare, either parent, or sister node spanning root
    	int randomNodeIdx = Randomizer.nextInt(tree.getNodeCount() - 1); // can't pick root
    	
    	System.out.println("randomNodeIdx=" + randomNodeIdx);
    	
        Node randomNode = tree.getNode(randomNodeIdx);
        Node parentNode = randomNode.getParent();
        Node nodeToCompare;
        if (parentNode.isRoot()) {
        	System.out.println("Spanning root when getting node to compare.");
        	Node left = parentNode.getChild(0);
        	Node right = parentNode.getChild(1);
        	int leftIdx = left.getNr();
        	int rightIdx = right.getNr();
        	
        	System.out.println("randomNodeIdx=" + randomNodeIdx + " leftIdx=" + leftIdx +  " rightIdx= " + rightIdx);
        	
        	if (leftIdx != randomNodeIdx) {
        		nodeToCompare = left;
        	} else { nodeToCompare = right; }
        }
        else {
        	nodeToCompare = parentNode;
        }
       
        int nodeToCompareIdx = nodeToCompare.getNr();
       
        System.out.println("Picked node to compare with idx: " + nodeToCompareIdx);
        
        // now to the colors
        int currentColor = colorAssignments[randomNodeIdx];
        int colorToCompare = colorAssignments[nodeToCompareIdx];
        int proposedColor;
        if (inactiveColorsSet.isEmpty()) {
        	proposedColor = currentColor; // merge because everyone has a different color
        }
        else {
        	Iterator<Integer> iter = inactiveColorsSet.iterator();
            proposedColor = iter.next(); // get first inactive color
        }
        
        System.out.println("Current color: " + nodeToCompareIdx);
        System.out.println("Color to compare: " + nodeToCompareIdx);
        System.out.println("Proposed color: " + nodeToCompareIdx);
        
        // now deciding whether to merge or split
        double hastingsRatio;
        boolean merged = false;
        double currentNColors = activeColorsSet.size();
        System.out.println("Current # of active colors=" + currentNColors);
        if (currentColor != colorToCompare) {
        	mergeAndUpdateAssignments(colorAssignments, currentColor, colorToCompare, proposedColor);
        	hastingsRatio = (2*nSpp - 2 - currentNColors + 1) / currentNColors;
        	merged = true;
        }
        else {
        	splitAndUpdateAssignments(colorAssignments, randomNode, currentColor, proposedColor);
        	hastingsRatio = (currentNColors + 1) / (2*nSpp - 2 - currentNColors);
        }
        
        // finally update nColors
        if (merged) { nColorInput.get().setValue((int) currentNColors - 1);
        } else { nColorInput.get().setValue((int) currentNColors + 1); }
        
        System.out.println("HR=" + hastingsRatio);
        System.out.println("After proposal: " + Arrays.toString(colorAssignmentInput.get().getValues()) + " nColors=" + nColorInput.get().getValue());
        return Math.log(hastingsRatio);
        //return Math.log(1.0);
    }
    
    private void mergeAndUpdateAssignments(Integer[] colorAssignments, int currentColor, int colorToCompare, int proposedColor) {
    	System.out.println("Merging.");
    	
    	// if we want to paint the root as well, change to .length-1
    	for (int i=0; i<(colorAssignments.length); ++i) {
    		if (colorAssignments[i] == currentColor || colorAssignments[i] == colorToCompare) {
    			colorAssignmentInput.get().setValue(i, proposedColor);
    		}
    	}
    }
    
    private void splitAndUpdateAssignments(Integer[] colorAssignments, Node randomNode, int currentColor, int proposedColor) {
    	System.out.println("Splitting.");
    	
    	// first, paint this node's subtending branch
    	int randomNodeIdx = randomNode.getNr();
    	colorAssignmentInput.get().setValue(randomNodeIdx, proposedColor);
    	
    	// now recur down the subtree
    	recursivelyPaintSubTree(colorAssignments, randomNode, currentColor, proposedColor);
    }
    
    private void recursivelyPaintSubTree(Integer[] colorAssignments, Node aNode, int currentColor, int proposedColor) {
    	int nodeIdx = aNode.getNr();
    	
    	// no matter what type of node, we paint if the colors match
    	if (colorAssignments[nodeIdx] == currentColor) {
    		colorAssignmentInput.get().setValue(nodeIdx, proposedColor);
    	}
    	
    	// if either leaf or internal node w/ different color, return
    	else if (aNode.isLeaf()) {
    		return;
    	}
    	
    	// if internal node still has matching color, keep recurring until either tip or mismatch between colors
    	else {
    		Node left = aNode.getChild(0);
        	Node right = aNode.getChild(1);
        	
        	recursivelyPaintSubTree(colorAssignments, left, currentColor, proposedColor);
        	recursivelyPaintSubTree(colorAssignments, right, currentColor, proposedColor);
    	}
    }
}