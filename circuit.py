import copy
import math
from typing import List
from pauli_path import PauliPath

class Circuit:

    def __init__(self, num_qubits:int, num_layers:int, l:int, gate_pos:List[List[tuple]]):
        """
        Initiate all the lists of layers for a given weight configuration.
        """
        self.num_qubits = num_qubits
        self.num_layers = num_layers
        self.gate_pos = gate_pos

        self.max_weight = l # upper bound on a Pauli path's Hamming weight

        self.weight_combos = []

        self.enumerate_weights([], self.max_weight-self.num_layers, self.num_layers)
        self.init_pauli_paths()
    

    def init_pauli_paths(self):
        self.pauli_paths = []
        for weight_combo in self.weight_combos:
            self.pauli_paths.append(PauliPath(self.num_qubits, weight_combo, self.gate_pos))


    """
    This function determines all solutions to the equation w_1 + w_2 + ... + w_d = k, where w_i >= 1,
    and we take into account the restrictions of the circuit architecture and legal Pauli path requirements

    This function takes three arguments and recursively constructs our list of all weight combinations.

    Args:
        weight_list (List[int]) : List of the weight combination we're currently working with
        wiggle_room (int) : Number of layers that can have weight greater than 1
        num_layers_left (int) : Number of layers waiting for weight

    Returns:
        void: The function appends to weight_combos, which affects the list in the calling function
    """
    def enumerate_weights(self, weight_list:List[int], wiggle_room:int, num_layers_left:int):

        # Base case: no more layers to add weight to
        if num_layers_left == 0:
            self.weight_combos.append(weight_list)
            return

        # Base case: no more wiggle room
        if wiggle_room == 0:
            if weight_list[-1] == 2: # Otherwise we would not meet the legal Pauli path requirements
                for i in range(num_layers_left):
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
                max_extra_weight = wiggle_room+1
            else:
                min_extra_weight = max(0,weight_list[-1]-len(self.gate_pos[len(weight_list)-1])-1) # -1 since we add 1 in the for loop
                max_extra_weight = min(len(self.gate_pos[len(weight_list)-1])*2, weight_list[-1]*2, wiggle_room+1)
            
            for i in range(min_extra_weight, max_extra_weight):
                weight_list_copy = copy.deepcopy(weight_list) # Copy of weight_list to avoid overwriting
                weight_list_copy.append(i+1) # +1 since every weight is guaranteed to be least 1
                self.enumerate_weights(weight_list_copy, wiggle_room-i, num_layers_left-1)
            return