import pdb
import copy
import math
from typing import List
from pauli_operator import PauliOperator
from pauli_path_traversal import PauliPathTraversal

class CircuitSim:

    def __init__(self, num_qubits:int, l:int, gate_pos:List[List[tuple]]):
        """
        Initiate all the lists of layers for a given weight configuration.
        """
        self.num_qubits = num_qubits
        self.num_op_layers = len(gate_pos)+1
        self.gate_pos = gate_pos

        self.max_weight = l # upper bound on a Pauli path's Hamming weight

        self.weight_combos = []

        self.enumerate_weights([], self.max_weight-self.num_op_layers, self.num_op_layers)
        self.init_pauli_paths()
        self.travs_to_list()
    
    def travs_to_list(self):
        self.pauli_path_list = []
        for pauli_path_trav in self.pauli_paths:
            self.pauli_path_list.append(self.trav_to_list(pauli_path_trav))

    def trav_to_list(self,trav:PauliPathTraversal):
        trav_list = []
        for sibs in trav.layers[0].forward_sibs.values():
            for pauli_op in sibs:
                self.pauli_op_hopping(trav_list, [], pauli_op)
        return trav_list

    def pauli_op_hopping(self,trav_list:List[List[List[str]]], partial_pauli_path:List[List[str]], pauli_op:PauliOperator):
        partial_pauli_path.append(pauli_op.operator)

        # Base case: Reached last Pauli operator layer of the circuit
        if pauli_op.next_ops == None:
            trav_list.append(partial_pauli_path)
            return
    
        for i in range(len(pauli_op.next_ops)):
            partial_pauli_path_copy = copy.deepcopy(partial_pauli_path)
            self.pauli_op_hopping(trav_list, partial_pauli_path_copy, pauli_op.next_ops[i])

    def init_pauli_paths(self):
        self.pauli_paths = []
        for weight_combo in self.weight_combos:
            self.pauli_paths.append(PauliPathTraversal(self.num_qubits, weight_combo, self.gate_pos))


    
    """
    This function determines all solutions to the equation w_1 + w_2 + ... + w_d = k, where w_i >= 1,
    and we take into account the restrictions of the circuit architecture and legal Pauli path requirements

    This function takes three arguments and recursively constructs our list of all weight combinations.

    Args:
        weight_list (List[int]) : List of the weight combination we're currently working with
        wiggle_room (int) : Number of layers that can have weight greater than 1
        num_op_layers_left (int) : Number of layers waiting for weight

    Returns:
        void: The function appends to weight_combos, which affects the list in the calling function
    """
    def enumerate_weights(self, weight_list:List[int], wiggle_room:int, num_op_layers_left:int):
        
        # Base case: no more layers to add weight to
        if num_op_layers_left == 0:
            self.weight_combos.append(weight_list)
            return

        # Base case: no more wiggle room
        if wiggle_room == 0:
            if weight_list[-1] <= 2: # Otherwise we would not meet the legal Pauli path requirements
                for i in range(num_op_layers_left):
                    weight_list.append(1) # All remaining layers get weight 1, since there's no wiggle room
                self.weight_combos.append(weight_list)
            return

        # recursive case: still have wiggle room and layers to add weight to
        else:
            # Lower bound on next weight: prior weight / 2, 
            # which is the case where all gates between prior and current layer have RR as input
            # upper bound on next weight: prior weight * 2, 
            # which is the case where all gates between prior and current layer have IR or RI as input
            # the range represents how much we are adding in addition to the 1 guranteed every layer
            # need new_weight <= 2^{n-1}-3 if our remaining weight <= 2^{n}-n-3
            if len(weight_list) == 0:
                min_extra_weight = 0
                max_extra_weight = min(self.num_qubits,wiggle_room+1) # +1 since this will be used as the excluded upper bound
            else:
                if len(self.gate_pos[len(weight_list)-1]) < math.floor(weight_list[-1]/2): # there aren't enough gates to halve the Hamming weight at the prior layer
                    min_extra_weight = weight_list[-1]-len(self.gate_pos[len(weight_list)-1])-1 # -1 since we add 1 in the for loop
                else:
                    min_extra_weight = math.ceil(weight_list[-1]/2)-1 # -1 since we always append +1 in our for loop 
                max_extra_weight = min(weight_list[-1]+len(self.gate_pos[len(weight_list)-1]), weight_list[-1]*2, wiggle_room+1,self.num_qubits)
            
            for i in range(min_extra_weight, max_extra_weight):
                weight_list_copy = copy.deepcopy(weight_list) # Copy of weight_list to avoid overwriting
                weight_list_copy.append(i+1) # +1 since every weight is guaranteed to be least 1
                self.enumerate_weights(weight_list_copy, wiggle_room-i, num_op_layers_left-1)
            return