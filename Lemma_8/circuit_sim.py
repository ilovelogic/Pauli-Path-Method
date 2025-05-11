import pdb
import copy
import math
from typing import List
from pauli_operator import PauliOperator
from pauli_path_trav import PauliPathTrav
from sibling_ops import SiblingOps

class CircuitSim:
    """
    Initiates a list of all PauliPathTrav objects that satisfy the constraints of the circuit and 
    upper bound on total Hamming weight of a Pauli path.
    Also fills out a list of all list representations of legal Pauli paths given the circuit architecture
    and Hamming weight upper bound.
    """
    def __init__(self, num_qubits:int, l:int, gate_pos:List[List[tuple]]):

        if not self.valid_gate_pos(num_qubits,gate_pos):
            print(gate_pos)
            raise ValueError
        
        self.num_qubits = num_qubits
        self.num_op_layers = len(gate_pos)+1
        self.gate_pos = gate_pos
        self.max_weight = l # upper bound on a Pauli path's Hamming weight

        if self.max_weight < self.num_op_layers: # we cannot make a valid Pauli path
            raise ValueError

        self.weight_combos = []
        # Fills out self.weight_combos with all weight configurations of legal Pauli paths
        # given the circuit and Hamming weight upper bound
        self.enumerate_weights([], self.max_weight-self.num_op_layers, self.num_op_layers)

        self.init_pauli_paths() # Adds the PauliPathTrav that matches each weight combo to self.pauli_path_travs
        self.travs_to_list()

        self.build_xyz_tree()
        return
        self.tree_to_lists()

    @staticmethod
    def valid_gate_pos(num_qubits:int, gate_pos:List[List[tuple]]):
        for gate_pos_layer in gate_pos:
            gate_pos_layer_set = set()
            for pos in gate_pos_layer:
                # Check that the positions reference qubits in the circuit
                if pos[0] < 0 or pos[0] > num_qubits-1:
                    return False
                if pos[1] < 0 or pos[1] > num_qubits-1:
                    return False
                # Check that there are no qubits that are input to two different gates at a gate layer
                if pos[0] in gate_pos_layer_set or pos[1] in gate_pos_layer_set:
                    return False
                gate_pos_layer_set.add(pos[0])
                gate_pos_layer_set.add(pos[1])
        return True
            


    """
    This function determines all solutions to the equation w_1 + w_2 + ... + w_d = k, where w_i >= 1,
    and we take into account the restrictions of the circuit architecture and legal Pauli path requirements

    This function takes three arguments and recursively constructs our list of all weight combinations.

    Args:
        weight_list (List[int]) : List of the weight combination we're currently working with
        wiggle_room (int) : Number of remaining Pauli operator layers that can have weight greater than 1
        num_op_layers_left (int) : Number of Pauli operator layers waiting to be assigned their weight

    Returns:
        void: The function appends each finished weight_list to the calling object's weight_combos
    """
    def enumerate_weights(self, weight_list:List[int], wiggle_room:int, num_op_layers_left:int):
        
        # Base case: no more layers to add weight to
        if num_op_layers_left == 0:
            self.weight_combos.append(weight_list)
            return

        # Base case: no more wiggle room
        if wiggle_room == 0:
            if len(weight_list) == 0 or weight_list[-1] <= 2: # Otherwise we would not meet the legal Pauli path requirements,
                # since we would decrease the Hamming weight by more than half between two layers
                for i in range(num_op_layers_left):
                    weight_list.append(1) # All remaining layers get weight 1, since there's no wiggle room
                self.weight_combos.append(weight_list)

            return

        # Recursive case: still have wiggle room and Pauli operator layers to add weight to
        else:
            if len(weight_list) == 0:
                min_extra_weight = 0
                max_extra_weight = min(self.num_qubits,wiggle_room+1) # +1 since this will be used as the excluded upper bound
            else:
                # Lower bound on next weight: The case where all the 'R's in the prior layer are grouped
                # to be 'RR' inputs to gates that have 'RI' or 'IR' as output
               
                # Checks if there aren't enough gates to halve the Hamming weight at the prior layer
                # or if the Hamming weight of the prior layer can be reduced by 1 by every gate
                if len(self.gate_pos[len(weight_list)-1]) < math.floor(weight_list[-1]/2): 
                    min_extra_weight = weight_list[-1]-len(self.gate_pos[len(weight_list)-1])-1 # -1 since we add 1 in the for loop
                else:
                    min_extra_weight = math.ceil(weight_list[-1]/2)-1 # -1 since we always append +1 in our for loop 
                # Upper bound on next weight: The case where all 'R's in the prior layer are sent with an 'I' into gates
                # that have the output of 'RR'
                # We can never more than double our prior weight, so if the prior weight is 1 and we have 10 gates
                # the best we can do is get 2
                # We don't subtract 1 like we did for min_extra_weight since max_extra_weight is the excluded upper bound
                max_extra_weight = min(weight_list[-1]+len(self.gate_pos[len(weight_list)-1]), weight_list[-1]*2, wiggle_room+1,self.num_qubits)
            
            for i in range(min_extra_weight, max_extra_weight):
                weight_list_copy = copy.deepcopy(weight_list) # Copy of weight_list to avoid overwriting
                weight_list_copy.append(i+1) # +1 since every weight is guaranteed to be least 1
                self.enumerate_weights(weight_list_copy, wiggle_room-i, num_op_layers_left-1)
            return

    # Initiates the list of all PauliPathTravs that fit the circuit architecture and upper bound on Hamming weight
    def init_pauli_paths(self):
        self.pauli_path_travs = []
        for weight_combo in self.weight_combos:
            self.pauli_path_travs.append(PauliPathTrav(self.num_qubits, weight_combo, self.gate_pos))
            

    
    # Translates PauliPathTrav objects into a list of paths, still in 'R', 'N', 'P', and 'I'
    def travs_to_list(self):
        self.rnp_pauli_paths = []
        for pauli_path_trav in self.pauli_path_travs:
            self.rnp_pauli_paths.append(self.trav_to_list(pauli_path_trav))

    def trav_to_list(self,trav:PauliPathTrav):
        trav_list = []
        for sibs in trav.layers[0].forward_sibs.values():
            for pauli_op in sibs:
                self.pauli_op_hopping(trav_list, [], pauli_op)
        return trav_list

    def pauli_op_hopping(self,trav_list:List[List[PauliOperator]], partial_pauli_path:List[PauliOperator], pauli_op:PauliOperator):
        partial_pauli_path.append(pauli_op)

        # Base case: Reached last Pauli operator layer of the circuit
        if pauli_op.next_ops == None:
            trav_list.append(partial_pauli_path)
            return
    
        for i in range(len(pauli_op.next_ops)):
            partial_pauli_path_copy = copy.deepcopy(partial_pauli_path)
            self.pauli_op_hopping(trav_list, partial_pauli_path_copy, pauli_op.next_ops[i])

    # Translates each Pauli path list in 'R', 'N', 'P', and 'I' to construct a tree 
    # structure, with its branching representing valid selections of 'X', 'Y', and 'Z'
    def build_xyz_tree(self):
        self.sib_op_heads = []
        for list_of_paths in self.rnp_pauli_paths:
            for path in list_of_paths:
                first_op_list = self.rn_to_z(path[0]) # returns a list with single element, 'I' 'Z' version of path[00]

                sib_op = SiblingOps(first_op_list, 1, path)
                self.sib_op_heads.append(sib_op)

    # Turns each tree into seperate lists representing Pauli paths
    def tree_to_lists(self):
        self.xyz_pauli_paths = []
        
        for sib_op_head in self.sib_op_heads:
            sib_op_paths = [[]]
            self.branch(sib_op_head,sib_op_paths)
    

    def branch(self, cur_sib:SiblingOps, pauli_paths_in_womb:List[List[PauliOperator]]):
        if (cur_sib.next_sibs == None):
            for pauli_path in pauli_paths_in_womb:
                pauli_path.append(cur_sib.pauli_ops[0]) # only the next pauli operator is the one with all 'I's and 'Z's
                self.xyz_pauli_paths.append(pauli_path) # adds all the completed pauli paths
        else:
            for pauli_op in cur_sib.pauli_ops:
                for next_sib in cur_sib.next_sibs:
                    branched_pauli_paths = copy.deepcopy(pauli_paths_in_womb)
                    for pauli_path in branched_pauli_paths:
                        pauli_path.append(pauli_op)
                    self.branch(next_sib,branched_pauli_paths)
            

    @staticmethod
    # The first Pauli operator in a Pauli path can only be a tensor of 'I's and 'Z's
    def rn_to_z(first_op:PauliOperator):
        first_op_list = [copy.deepcopy(first_op)]

        for i in range(len(first_op_list[0].operator)):
            if first_op_list[0].operator[i] == 'R' or first_op_list[0].operator[i] == 'N':
                first_op_list[0].operator[i] = 'Z'
            if first_op_list[0].operator[i] == 'P':
                return [] # we shouldn't encounter a 'P' in the first layer
            
        return first_op_list # valid operator possible