from __future__ import annotations
from typing import List
import copy

class PauliOperator:

<<<<<<< HEAD
    def __init__(self, operator:List[str], backward_ops:List[PauliOperator] = None, forward_ops:List[PauliOperator] = None):
=======
    def __init__(self, operator:List[str], prev_ops:List[PauliOperator] = None, next_ops:List[PauliOperator] = None):
>>>>>>> cb96eb3 (updated former class PathLayer to be PauliOperator. rewrote weight_to_operators function to take into account the weight of non-gate qubits.)
        """
        Initialize.
        """
                  
        self.operator = operator
<<<<<<< HEAD
        self.backward_ops = backward_ops
        self.forward_ops = forward_ops

    def weight_to_operators(self, next_weight:int, pos_to_fill:List[tuple], backward:int):
=======
        self.prev_ops = prev_ops
        self.next_ops = next_ops

    def weight_to_layers(self, next_weight:int, pos_to_fill:List[tuple], backward:int):
>>>>>>> cb96eb3 (updated former class PathLayer to be PauliOperator. rewrote weight_to_operators function to take into account the weight of non-gate qubits.)
        next_gate_weight = next_weight

        unordered_pos_to_fill = {pos for gate_pos in pos_to_fill for pos in gate_pos}
        for i in range(len(self.operator)):
            if i not in unordered_pos_to_fill and self.operator[i] == 'R':
                next_gate_weight -= 1 # every non-gate qubit that is non-identity takes from our 
                # overall Hamming weight available to gate qubits


        self.list_alloc = self.list_allocs(len(pos_to_fill),next_gate_weight)
        
        sibs = []

        for i in range(self.list_alloc[len(pos_to_fill)][next_gate_weight]): 
            sibs.append(PauliOperator(copy.deepcopy(self.operator))) # copies layer_str and uses to initialize PauliOperators

        num_RRs = next_gate_weight - len(pos_to_fill) # number of RRs we can use to fill in the layer

<<<<<<< HEAD
        # no way to make a valid layer, given the weights
        if num_RRs < 0 or num_RRs > len(pos_to_fill):
            print("Can't make it work, baby!")
            
        else:
            if (backward):
                self.backward_ops = sibs
                self.add_gate_input(self.backward_ops, num_RRs, pos_to_fill, 0)
            else:
                self.forward_ops = sibs
                self.add_gate_input(self.forward_ops, num_RRs, pos_to_fill, 0)
=======
        # distingush between rrs and irs
        if num_RRs < 0 or num_RRs > len(pos_to_fill):
            print("Could not satisfy user input with a valid layer")
            # no way to make a valid layer, given the weights
        else:
            if (backward):
                self.prev_ops = sibs
                self.add_gate_input(self.prev_ops, num_RRs, pos_to_fill, 0)
            else:
                self.next_ops = sibs
                self.add_gate_input(self.next_ops, num_RRs, pos_to_fill, 0)
>>>>>>> cb96eb3 (updated former class PathLayer to be PauliOperator. rewrote weight_to_operators function to take into account the weight of non-gate qubits.)


    """
    This function determines the number of entries we need to allocate in our list 
    for our recursive call with a given number of non-identity gate positions and the Hamming weight.

    It takes two arguments and uses a dynamic programming approach to fill out a 2D array.

    Args:
        num_p (int) : the number of non-identity gate positions we have to fill
        num_w (int) : the Hamming weight we have to spread across the gate positions
        
    Returns:
        List[List[int]] : a 2D array where the (i,j) entry stores the number of ways we
        can fill i gates with non-identity I/O using a total Hamming weight of j
    """
        
    @staticmethod
    def list_allocs(num_p:int, num_w:int):
        list_alloc = [[0 for _ in range(num_w+1)] for _ in range(num_p+1)]

        # base cases
        for w in range(0, num_w):
            list_alloc[0][w] = 0 # no gates, no way to have a layer

        for p in range(0, num_p): # no weight, only one way to have a layer
            list_alloc[p][0] = 1
        
        for i in range(1, num_p+1):
            for j in range(1,num_w+1):   #     IR                     RI                     RR
                if (j > i):
                    list_alloc[i][j] = list_alloc[i-1][j-1] + list_alloc[i-1][j-1] + list_alloc[i-1][j-2]
                else:
                    list_alloc[i][j] = list_alloc[i-1][j-1] + list_alloc[i-1][j-1]
            
        return list_alloc



    """
    This function determines all possible layers at a depth one away from this layer,
    given this layer and the weight of the new layer.

    It takes five arguments and uses helper function add_gate_input to build a list of all gate output/input combos

    Args:
        prior_layer (str) : the sting representing the prior layer of the circuit,
        with each character representing a qubit
        this_weight (int) : the Hamming weight of the new layer
        n (int) : the number of qubits
        pos_list (List[tuple]) : each tuple contains the indices of the two inputs to one of the gates


    Returns:
        void : the function updates self.layers to be the list of all valid layers at this depth in the circuit, 
    with each layer represented as a list of strings and each string being 'I', 'X', 'Y', or 'X'.
    """

    @staticmethod
    def append_to_layers(sibs:List[PauliOperator],indices:tuple, strs:tuple, r_start:int, r_end:int):
        ind1, ind2 = indices
        str1, str2 = strs
        for i in range(r_start, r_end):
            sibs[i].operator[ind1] = str1
            sibs[i].operator[ind2] = str2


    def add_gate_input(self, sibs:List[PauliOperator], num_RRs:int, pos_to_fill:List[tuple], r_start:int):

        if (len(pos_to_fill) == 0):
            return
        
        if (len(pos_to_fill) == num_RRs): # no more wiggle room, we must fill all remaining gate inputs with RR     
            cur_pos = pos_to_fill.pop(0)
            rr_start = r_start
            rr_end = r_start+self.list_alloc[len(pos_to_fill)+1][len(pos_to_fill)+num_RRs+1]
            self.append_to_layers(sibs, cur_pos, ('R','R'), rr_start, rr_end) # copy of layers with 'RR' added to all strs
            self.add_gate_input(sibs, num_RRs-1, list(pos_to_fill), rr_start) # one less RR to use,
            # next call will handle adding the next RR
            return

        else:
            cur_pos = pos_to_fill.pop(0)

            ir_start = r_start
            ir_end = r_start + self.list_alloc[len(pos_to_fill)][len(pos_to_fill)+num_RRs]
            self.append_to_layers(sibs, cur_pos, ('I','R'), ir_start, ir_end) # copy of layers with 'IR' added to all strs

            ri_start = r_start + self.list_alloc[len(pos_to_fill)][len(pos_to_fill)+num_RRs]
            ri_end = r_start + 2*self.list_alloc[len(pos_to_fill)][len(pos_to_fill)+num_RRs]
            self.append_to_layers(sibs, cur_pos, ('R','I'), ri_start, ri_end) # copy of layers with 'RI' added to all strs

            if (num_RRs != 0):
                rr_start = r_start + 2*self.list_alloc[len(pos_to_fill)][len(pos_to_fill)+num_RRs]
                rr_end = r_start+self.list_alloc[len(pos_to_fill)+1][len(pos_to_fill)+num_RRs+1]
                self.append_to_layers(sibs, cur_pos, ('R','R'), rr_start, rr_end) # copy of layers with 'RR' added to all strs
                self.add_gate_input(sibs, num_RRs-1, list(pos_to_fill), rr_start) 

            self.add_gate_input(sibs, num_RRs, list(pos_to_fill), ir_start) 
            self.add_gate_input(sibs, num_RRs, list(pos_to_fill), ri_start) 
            
            return