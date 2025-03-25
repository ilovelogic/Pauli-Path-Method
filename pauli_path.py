import pdb # for debugging
import copy
from collections import defaultdict
from typing import List, DefaultDict
from itertools import combinations
from pauli_operator import PauliOperator
from layer import Layer

class PauliPath:

    """
    Initiate all the lists of layers for a given weight configuration.
    """
        
    def __init__(self, num_qubits:int, weight_combo:List[int],gate_pos:List[List[tuple]]):
        self.num_qubits = num_qubits
        self.depth = len(weight_combo) # weight_combo stores exactly one weight for each layer
        self.weight_combo = weight_combo
        self.gate_pos = gate_pos
        self.all_ops = [None]*self.depth # List of DefaultDicts of all operators possible at each layer
        
        if (num_qubits == 0 or self.depth == 0): # No circuit -> nothing to do
            return
        if (self.depth == 1): # Single layer -> generate for that layer, no propagation
            self.all_ops[0] = defaultdict(list)
            self.all_ops[0][(tuple([]), tuple(""))] =  self.unsorted_min_layer_ops(weight_combo[0])
            return

        min_layer_ops, pos_to_fill_b, pos_to_fill_f, min_depth = self.build_min_configs()
        self.all_ops[min_depth] = min_layer_ops

        if (min_depth-1 >= 0):
            min_next_sibs_b = self.propagate_next(min_layer_ops.backward_sibs, pos_to_fill_b, 1, min_depth)
            if (min_depth-2 >= 0):
                self.all_ops[min_depth-1] = Layer(gate_pos[min_depth-2], 1, min_next_sibs_b)
            else:
                self.all_ops[min_depth-1] = Layer()
                self.all_ops[min_depth-1].forward_sibs = min_next_sibs_b
        
        if (min_depth < self.depth-1):
            min_next_sibs_f = self.propagate_next(self.all_ops[min_depth].forward_sibs, pos_to_fill_f, 0, min_depth)
            if (min_depth+1 < self.depth-1):
                self.all_ops[min_depth+1] = Layer(gate_pos[min_depth+1], 0, min_next_sibs_f)
            else:
                self.all_ops[min_depth+1] = Layer()
                self.all_ops[min_depth+1].backward_sibs = min_next_sibs_f

        # Propagating backward
        for i in range(min_depth-2, -1, -1): # Goes from min_depth -1 to 0
            # Determines all the sibs that the prior depth can propagate backwards to 
            next_sibs_b = self.propagate_next(self.all_ops[i+1].backward_sibs, self.all_ops[i+1].pos_to_fill, 1, i+1)
            self.all_ops[i] = Layer(gate_pos[i-1], 1, next_sibs_b)

        # Propagating forward
        for i in range(min_depth+2, self.depth): # Goes from min_depth + 2 to self.depth-1
            next_sibs_f = self.propagate_next(self.all_ops[i-1].forward_sibs, self.all_ops[i-1].pos_to_fill, 0, i)
            self.all_ops[i] = Layer(gate_pos[min_depth+1], 0, next_sibs_f)

        return
    

    def build_min_configs(self):

        min_weight = float('inf') # Will get overwritten by first value we check

        for i in range(self.depth): # Iterates over all depths to find the one with lightest weight
            if (self.weight_combo[i] < min_weight):
                min_depth = i # Saves index of layer with min weight (so far)
                min_weight = self.weight_combo[i] # Saves min weight (so far)
                
        min_layers = self.unsorted_min_layer_ops(min_weight)

        if (min_depth == 0):
            min_layer_ops = Layer(self.gate_pos[min_depth], 0)
            min_layer_ops.check_qubits(min_layers)
            min_layer_ops.find_sibs(min_layers)
            # Since we are at the back-most layer, we return an empty list for the backward position list
            return min_layer_ops, [], min_layer_ops.pos_to_fill, min_depth
        
        min_layer_ops, pos_to_fill_b = self.min_backward(min_layers, min_depth)

        if (min_depth == self.depth-1): # Minimum depth is at the last layer of the circuit
            # since we are at the forward-most layer, we return an empty list for the forward position list
            return min_layer_ops, pos_to_fill_b, [], min_depth

        pos_to_fill_f = self.min_forward(min_layers, min_layer_ops, min_depth)

        return min_layer_ops, pos_to_fill_b, pos_to_fill_f, min_depth
    
    # Creates a list of all possible lists of indicies in increasing order 
    # of length min_weight using values from 0 to num_qubits-1
    # Ex. num_qubits = 3, min_weight = 2: [[0,1], [0,2], [1,2]]
    def unsorted_min_layer_ops(self, min_weight):
        arrangements = list(combinations(range(self.num_qubits), min_weight)) 
        min_layers = []
        for arrangement in arrangements: 
            temp = ['I'] * self.num_qubits # Generates list ['I','I',...] that has size num_qubits
            for index in arrangement:
                temp[index] = 'R' # Replaces all the indices in one of our arrangements with 'R'
                # Results in a layer with Hamming weight = min_weight
            min_layers.append(PauliOperator(temp))
        return min_layers

    def min_backward(self,min_layers,min_depth):
        # Setting min_layer_ops' backward_sibs,
        # which is the list of lists of sibling PauliOperators,
        # where each individual list is a grouping of PauliOperators at depth min_depth 
        # that propagate *backward* to the same list of PauliOperators
        min_layer_ops = Layer(self.gate_pos[min_depth-1], 1)
        min_layer_ops.check_qubits(min_layers)
        min_layer_ops.find_sibs(min_layers)
        pos_to_fill_b = min_layer_ops.pos_to_fill # For all Layers except for the one at min_depth
        # we are only propagating in one direction,
        # so we usually do not need to save the Layer's pos_to_fill in a seperate var
        # However, for min_depth, we are propagating forward and backward,
        # so we need to save the backward pos_to_fill so we do not lose it
        # when we use min_layer_ops to propagate forward and pos_to_fill is set to the forward version
        return min_layer_ops, pos_to_fill_b

    def min_forward(self,min_layers,min_layer_ops,min_depth): 
        # Setting min_layer_ops' forward_sibs,
        # which is the list of lists of sibling PauliOperators,
        # where each individual list is a grouping of PauliOperators at depth min_depth 
        # that propagate *forward* to the same list of PauliOperators
        min_layer_ops.gate_pos = self.gate_pos[min_depth] # for setting forward_sibs
        min_layer_ops.backward = 0
        min_layer_ops.check_qubits(min_layers)
        min_layer_ops.find_sibs(min_layers)
        pos_to_fill_f = min_layer_ops.pos_to_fill # We technically do not need a seperate var since we could
        # just reference min_configs.pos_to_fill
        # The reason for the new var is simply for readability in the init function
        return pos_to_fill_f

    """
    This function takes in a list of forward or backward sibling operators at a Layer,
    and determines the new sibling operators that each sibling operators list in the input list propagate to

    It takes four arguments and uses helper function weight_to_operaters from the PauliOperator class
    to get the sibling operators that all the PauliOperators in any given sibling operators of the input
    propagate to

    Args:
        all_sibs (DefaultDict[tuple, List[PauliOperator]]) : the hash map of all sibling operators at the prior depth of the circuit
        pos_to_fill (DefaultDict[PauliOperator,List]) : a hash map of the gate positions with non-identity I/O between the prior depth and next depth
        backward (int) : 1 if we are propagating backward, 0 if we are propagating forward
        depth (int) : the new depth for which we are getting sibling operators

    Returns:
        DefaultDict[tuple,List[PauliOperator]]: a list of all the sibling operators that the prior depth's 
        sibling operators propagate to
    """
    def propagate_next(self, all_sibs:DefaultDict[tuple, List[PauliOperator]], pos_to_fill:DefaultDict[PauliOperator,List], backward:int, depth:int):
        new_sib_ops = defaultdict(list)

        for sib_ops in all_sibs.values(): # Iterates over each sibling operators list at the prior depth 
            # All PauliOperators in a sibling operators have the same sibling operators they propagate to,
            # so we only need to get the next sibling operators once per each sibling operators in our list
            sib_ops[0].weight_to_operators(self.weight_combo[depth-backward], copy.deepcopy(pos_to_fill[sib_ops[0]]), backward)
        if (backward):
            for identifier, sib_ops in all_sibs.items():
                for op in sib_ops:
                    # Setting all PauliOperators in each sibling operators
                    # to have the same sibling operators that they propagate to 
                    op.backward_ops = sib_ops[0].backward_ops
                for backward_op in sib_ops[0].backward_ops: # Allowing for traversal
                    backward_op.forward_ops = sib_ops # in the reverse direction
                if (sib_ops[0].backward_ops != []): # Checks if propagation was successful
                    new_sib_ops[identifier] = sib_ops[0].backward_ops
        else:
            for identifier, sib_ops in all_sibs.items():
                for op in sib_ops:
                    op.forward_ops = sib_ops[0].forward_ops
                for forward_op in sib_ops[0].forward_ops:
                    forward_op.backward_ops = sib_ops
                if (sib_ops[0].forward_ops != []):
                    new_sib_ops[identifier] = sib_ops[0].forward_ops

        return new_sib_ops