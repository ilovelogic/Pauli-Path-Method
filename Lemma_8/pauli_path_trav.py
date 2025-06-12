import pdb
import copy
from collections import defaultdict
from typing import List, DefaultDict
from itertools import combinations
from pauli_operator import PauliOperator
from pauli_op_layer import PauliOpLayer

class PauliPathTrav:

    """
    Initiate all the lists of layers for a given weight configuration.
    """
        
    def __init__(self, num_qubits:int, weight_combo:List[int],gate_pos:List[List[tuple]]):

        self.num_qubits = num_qubits
        self.num_op_layers = len(weight_combo) # weight_combo stores exactly one weight for each layer
        self.weight_combo = weight_combo
        self.gate_pos = gate_pos
        self.layers = [None]*self.num_op_layers # List of PauliOpLayers of all operators possible at each layer
        
        if (num_qubits == 0 or self.num_op_layers == 0): # No circuit -> nothing to do
            return
        if (self.num_op_layers == 1): # Single layer -> generate for that layer, no propagation
            self.layers[0] = defaultdict(list)
            self.layers[0][(tuple([]), tuple(""))] =  self.unsorted_min_layer_ops(weight_combo[0])
            return

        min_layer_ops, pos_to_fill_b, pos_to_fill_f, min_index = self.build_min_configs()
        self.layers[min_index] = min_layer_ops

        if (min_index-1 >= 0): # min layer has a prior layer
            min_prior_sibs = self.propagate_next(min_layer_ops.backward_sibs, pos_to_fill_b, 1, min_index-1)
            if (min_index-2 >= 0): # min layer's prior layer has a prior layer
                self.layers[min_index-1] = PauliOpLayer(gate_pos[min_index-2], 1, min_prior_sibs)
            else:
                self.layers[min_index-1] = PauliOpLayer()
                self.layers[min_index-1].forward_sibs = min_prior_sibs

        if (min_index < self.num_op_layers-1):
            min_next_sibs_f = self.propagate_next(self.layers[min_index].forward_sibs, pos_to_fill_f, 0, min_index+1)
            if (min_index+1 < self.num_op_layers-1):
                self.layers[min_index+1] = PauliOpLayer(gate_pos[min_index+1], 0, min_next_sibs_f)
            else:
                self.layers[min_index+1] = PauliOpLayer()
                self.layers[min_index+1].backward_sibs = min_next_sibs_f

        # Propagating backward
        for i in range(min_index-2, -1, -1): # Goes from min_index -1 to 0
            # Determines all the sibs that the prior layer can propagate backwards to 
            
            next_sibs_b = self.propagate_next(self.layers[i+1].backward_sibs, self.layers[i+1].pos_to_fill, 1, i)
            self.layers[i] = PauliOpLayer(gate_pos[i-1], 1, next_sibs_b)

        # Propagating forward
        for i in range(min_index+2, self.num_op_layers): # Goes from min_index + 2 to self.num_op_layers-1
            next_sibs_f = self.propagate_next(self.layers[i-1].forward_sibs, self.layers[i-1].pos_to_fill, 0, i)
            self.layers[i] = PauliOpLayer(gate_pos[i-1], 0, next_sibs_f)

        return
    
    """
    This function finds the PauliOpLayer of the circuit with the lowest Hamming weight and
    determines all possible PauliOperators at that layer. It also seperately stores the 
    non-identity  gate positions between the min weight layer and its prior layer and 
    the min layer and its next layer.

    It uses the helper function unsorted_min_layer_ops to generate all possibile PauliOperators at the min layer.

    Args:
        self (PauliPathTrav) : The PauliPathTrav for which we need to find and initiate the min layer

    Returns:
        List[PauliOperator]: min_layer_ops, 
        List[tuple]: pos_to_fill_b, 
        List[tuple]: pos_to_fill_f, 
        int: min_index
    """

    def build_min_configs(self):

        min_weight = float('inf') # Will get overwritten by first value we check

        for i in range(self.num_op_layers): # Iterates over all layers to find the one with lightest weight
            if (self.weight_combo[i] < min_weight):
                min_index = i # Saves index of layer with min weight (so far)
                min_weight = self.weight_combo[i] # Saves min weight (so far)
                
        min_layer_ops_list = self.unsorted_min_layer_ops(min_weight)

        if (min_index == 0):
            min_layer_ops = PauliOpLayer(self.gate_pos[min_index], 0)
            min_layer_ops.check_qubits(min_layer_ops_list)
            min_layer_ops.find_sibs(min_layer_ops_list)
            # Since we are at the back-most layer, we return an empty list for the backward position list
            return min_layer_ops, [], min_layer_ops.pos_to_fill, min_index
        
        min_layer_ops, pos_to_fill_b = self.min_backward(min_layer_ops_list, min_index)

        if (min_index == self.num_op_layers-1): # Minimum layer is at the last layer of the circuit
            # since we are at the forward-most layer, we return an empty list for the forward position list
            return min_layer_ops, pos_to_fill_b, [], min_index

        pos_to_fill_f = self.min_forward(min_layer_ops_list, min_layer_ops, min_index)

        return min_layer_ops, pos_to_fill_b, pos_to_fill_f, min_index
    
    def unsorted_min_layer_ops(self, min_weight:int):
        arrangements = list(combinations(range(self.num_qubits), min_weight)) # Creates a list of all possible 
        # lists of indicies in increasing order of length min_weight using values from 0 to num_qubits-1
        # Ex. num_qubits = 3, min_weight = 2: [[0,1], [0,2], [1,2]]
        min_layer_ops_list = []
        for arrangement in arrangements: 
            temp = ['I'] * self.num_qubits # Generates list ['I','I',...] that has size num_qubits
            for index in arrangement:
                temp[index] = 'R' # Replaces all the indices in one of our arrangements with 'R'
                # Results in a layer with Hamming weight = min_weight
            min_layer_ops_list.append(PauliOperator(temp))
        return min_layer_ops_list

    def min_backward(self,min_layer_ops_list:List[PauliOperator],min_index:int):
        # Setting min_layer_ops' backward_sibs,
        # which is the list of lists of sibling PauliOperators,
        # where each individual list is a grouping of PauliOperators at index min_index in the Pauli path
        # that propagate *backward* to the same list of PauliOperators
        min_layer_ops = PauliOpLayer(self.gate_pos[min_index-1], 1)
        min_layer_ops.check_qubits(min_layer_ops_list)
        min_layer_ops.find_sibs(min_layer_ops_list)
        pos_to_fill_b = min_layer_ops.pos_to_fill # For all Layers except for the one at min_index
        # we are only propagating in one direction,
        # so we usually do not need to save the Layer's pos_to_fill in a seperate var
        # However, for min_index, we are propagating forward and backward,
        # so we need to save the backward pos_to_fill so we do not lose it
        # when we use min_layer_ops to propagate forward and pos_to_fill is set to the forward version
        return min_layer_ops, pos_to_fill_b

    def min_forward(self,min_layer_ops_list:List[PauliOperator],min_layer_ops:PauliOpLayer,min_index:int): 
        # Setting min_layer_ops' forward_sibs,
        # which is the list of lists of sibling PauliOperators,
        # where each individual list is a grouping of PauliOperators at index min_index in the Pauli path 
        # that propagate *forward* to the same list of PauliOperators
        min_layer_ops.gate_pos = self.gate_pos[min_index] # for setting forward_sibs
        min_layer_ops.backward = 0
        min_layer_ops.check_qubits(min_layer_ops_list)
        min_layer_ops.find_sibs(min_layer_ops_list)
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
        all_sibs (DefaultDict[tuple, List[PauliOperator]]) : the hash map of all sibling operators at a layer of the Pauli path
        pos_to_fill (DefaultDict[PauliOperator,List]) : a hash map of the gate positions with non-identity I/O between the 
        PauliOpLayer associated with all_sibs and the PauliOpLayer we propagate to
        backward (int) : 1 if we are propagating backward, 0 if we are propagating forward
        index (int) : the new index in the Pauli path for which we are determining sibling operators

    Returns:
        DefaultDict[tuple,List[PauliOperator]]: a list of all the sibling operators that 
        the input PauliOpLayer's sibling operators propagate to
    """
    def propagate_next(self, all_sibs:DefaultDict[tuple, List[PauliOperator]], pos_to_fill:DefaultDict[PauliOperator,List], backward:int, index:int):
        new_sib_ops = defaultdict(list)

        for sib_ops in all_sibs.values(): # Iterates over each sibling operators list at the former index 
            # All PauliOperators in a sibling operators have the same sibling operators they propagate to,
            # so we only need to get the next sibling operators once per each sibling operators in our list
                sib_ops[0].weight_to_operators(sib_ops, self.weight_combo[index], copy.deepcopy(pos_to_fill[sib_ops[0]]), backward)
            
        if (backward):
            for identifier, sib_ops in all_sibs.items():
                for op in sib_ops:
                    # Setting all PauliOperators in each sibling operators
                    # to have the same sibling operators that they propagate to 
                    op.prior_ops = sib_ops[0].prior_ops
                for prior_op in sib_ops[0].prior_ops: # Allowing for traversal
                    prior_op.next_ops = sib_ops # in the reverse direction
                if (sib_ops[0].prior_ops != []): # Checks if propagation was successful
                    new_sib_ops[identifier] = sib_ops[0].prior_ops

        else:
            for identifier, sib_ops in all_sibs.items():
                for op in sib_ops:
                    op.next_ops = sib_ops[0].next_ops
                for next_op in sib_ops[0].next_ops:
                    next_op.prior_ops = sib_ops
                if (sib_ops[0].next_ops != []):
                    new_sib_ops[identifier] = sib_ops[0].next_ops

        return new_sib_ops