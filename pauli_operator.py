from __future__ import annotations
import pdb # for debuggingv

from typing import List
import copy

class PauliOperator:

    def __init__(self, operator:List[str], backward_ops:List[PauliOperator] = None, forward_ops:List[PauliOperator] = None):
        """
        Initialize.
        """
                  
        self.operator = operator
        self.backward_ops = backward_ops
        self.forward_ops = forward_ops

    """
    This function determines all possible operators at the Layer one depth away from this operator,
    that this operator can propagate to.

    It takes five arguments and uses helper function find_next_operators to build the propagation list.

    Args:
        self (PauliOperator) : The PauliOperator from which we are propagating
        next_weight (int) : The total Hamming weight of the next Layer
        pos_to_fill (List[tuple]) : The non-identity gate positions between the Layers
        backward (int) : 1 if we are propagating backward, 0 if we are propagating forward

    Returns:
        void : updates self.forward_ops or self.backward_ops (depending on the direction of propagation) 
        to be the list of all operators at the next layer of the circuit that this operator
        can propagate to.
    """
    def weight_to_operators(self, next_weight:int, pos_to_fill:List[tuple], backward:int):
        next_gate_weight = next_weight

        unordered_pos_to_fill = {pos for gate_pos in pos_to_fill for pos in gate_pos}
        for i in range(len(self.operator)):
            if i not in unordered_pos_to_fill and self.operator[i] == 'R':
                next_gate_weight -= 1 # Every non-gate qubit that is non-identity takes from our 
                # overall Hamming weight available to gate qubits

        num_RRs = next_gate_weight - len(pos_to_fill) # Number of RRs we can use to fill in the layer

        # No way to make a valid layer, given the weights
        if num_RRs < 0 or num_RRs > len(pos_to_fill):
            if (backward):
                self.backward_ops = []
            else:
                self.forward_ops = []
            return
        
        self.list_alloc = self.list_allocs(len(pos_to_fill),next_gate_weight)

        sibs = []
        for i in range(self.list_alloc[len(pos_to_fill)][next_gate_weight]): 
            sibs.append(PauliOperator(copy.deepcopy(self.operator))) # Copies layer_str and uses to initialize PauliOperators
    
        if (backward):
            self.backward_ops = sibs
            self.find_next_operators(self.backward_ops, num_RRs, pos_to_fill, 0)
        else:
            self.forward_ops = sibs
            self.find_next_operators(self.forward_ops, num_RRs, pos_to_fill, 0)

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
        if (num_p > num_w):
            return [] # No way to have num_p gates with non-identity I/O unless
            # we have at least num_p Hamming weight to spread across the gates

        list_alloc = [[0 for _ in range(num_w+1)] for _ in range(num_p+1)]

        # Base cases
        for w in range(1, num_w+1):
            list_alloc[0][w] = 0 # No positions to fill but positive weight -> no way to have a layer
        list_alloc[0][0] = 1 
        
        for i in range(1, num_p+1):
            for j in range(1,num_w+1):   #     IR                     RI                     RR
                if (j > i):
                    list_alloc[i][j] = list_alloc[i-1][j-1] + list_alloc[i-1][j-1] + list_alloc[i-1][j-2]
                else: # Can use only one non-identity Pauli per gate
                    list_alloc[i][j] = list_alloc[i-1][j-1] + list_alloc[i-1][j-1]
            
        return list_alloc

    @staticmethod
    def edit_operator(sibs:List[PauliOperator],indices:tuple, strs:tuple, r_start:int, r_end:int):
        ind1, ind2 = indices
        str1, str2 = strs
        for i in range(r_start, r_end):
            sibs[i].operator[ind1] = str1
            sibs[i].operator[ind2] = str2


    def find_next_operators(self, sibs:List[PauliOperator], num_RRs:int, pos_to_fill:List[tuple], r_start:int):

        if (len(pos_to_fill) == 0):
            return
        
        if (len(pos_to_fill) == num_RRs): # No more wiggle room, we must fill all remaining gate inputs with RR     
            cur_pos = pos_to_fill.pop(0)
            rr_start = r_start
            self.edit_operator(sibs, cur_pos, ('R','R'), rr_start, rr_start+1) # Copy of layers with 'RR' added to all operators
            self.find_next_operators(sibs, num_RRs-1, list(pos_to_fill), rr_start) # One less RR to use
            # Next call will handle adding the next RR
            return

        else:
            cur_pos = pos_to_fill.pop(0)

            ir_start = r_start
            ir_end = r_start + self.list_alloc[len(pos_to_fill)][len(pos_to_fill)+num_RRs]
            self.edit_operator(sibs, cur_pos, ('I','R'), ir_start, ir_end) # Copy of layers with 'IR' added to all operators

            ri_start = r_start + self.list_alloc[len(pos_to_fill)][len(pos_to_fill)+num_RRs]
            ri_end = r_start + 2*self.list_alloc[len(pos_to_fill)][len(pos_to_fill)+num_RRs]
            self.edit_operator(sibs, cur_pos, ('R','I'), ri_start, ri_end) # Copy of layers with 'RI' added to all operators

            if (num_RRs != 0):
                rr_start = r_start + 2*self.list_alloc[len(pos_to_fill)][len(pos_to_fill)+num_RRs]
                rr_end = r_start+self.list_alloc[len(pos_to_fill)+1][len(pos_to_fill)+num_RRs+1]
                self.edit_operator(sibs, cur_pos, ('R','R'), rr_start, rr_end) # Copy of layers with 'RR' added to all operators
                self.find_next_operators(sibs, num_RRs-1, list(pos_to_fill), rr_start) 

            self.find_next_operators(sibs, num_RRs, list(pos_to_fill), ir_start) 
            self.find_next_operators(sibs, num_RRs, list(pos_to_fill), ri_start) 
            
            return
        
    def r_to_xyz(self):
        self.r_pos = []
        for i in range(len(self.operator)):
            if self.operator[i] == 'R':
                self.r_pos.append(i)
        
        self.xyz_paulis = []
        for j in range(3 ** len(self.r_pos)):
            self.xyz_paulis.append(copy.deepcopy(self.operator))
        self.fill_in_r_pos('X',0,0)
        self.fill_in_r_pos('Y',0,3 ** (len(self.r_pos)-1))
        self.fill_in_r_pos('Z',0,2 * (3 ** (len(self.r_pos)-1)))

    def fill_in_r_pos(self, pauli:str, r_pos_i:int, start:int):
        if r_pos_i == len(self.r_pos):
            return
        for i in range(start, start + (3 ** (len(self.r_pos)-r_pos_i-1))):
            
            self.xyz_paulis[i][self.r_pos[r_pos_i]] = pauli
        self.fill_in_r_pos('X', r_pos_i+1,start)
        self.fill_in_r_pos('Y', r_pos_i+1,start+(3 ** (len(self.r_pos)-r_pos_i-2)))
        self.fill_in_r_pos('Z', r_pos_i+1,start+(2 * (3 ** (len(self.r_pos)-r_pos_i-2))))
