from collections import defaultdict
import copy
from typing import DefaultDict, List
from Path_Generation.pauli_operator import PauliOperator

class PauliOpLayer:

    def __init__(self, gate_pos:List[tuple]=None, backward:int=-1,pauli_ops:DefaultDict[tuple, List[PauliOperator]]=None):
        """
        Initiates the forward and backward lists of operators
        for propagating forward and backward at a given depth of the circuit
        """
        
        self.backward = backward # If backward = 1, propagating backward
        # Otherwise, propagating forward

        self.gate_pos = gate_pos # Gate positions between this Layer and:
        # If (backward): the previous Layer
        # Else: the next Layer

        if (pauli_ops is not None and self.gate_pos is not None and self.backward != -1):
            unsorted_pauli_ops = []

            # Builds list of all PauliOperators at this depth of the circuit
            # Turns DefaultDict into a 1D list
            for sibs in pauli_ops.values():
                for op in sibs:
                    unsorted_pauli_ops.append(op)

            self.check_qubits(unsorted_pauli_ops)

            if (self.backward):
                self.forward_rnp_sibs = pauli_ops # The input list in this case comes sorted according to
                # the lists that the Layer at the next depth propagated to when it propagated backward
            else: 
                self.backward_rnp_sibs = pauli_ops # The input list in this case comes sorted according to
                # the lists that the Layer at the previous depth propagated to when it propagated forward
            
            self.group_sibs(unsorted_pauli_ops) # Need to sort this list according to which ops at this depth propagate together 
            # either backward (if self.backward) or forward (else)

    """
    This function determines which gate positions of layers have non-identity I/O,
    and keeps track of each layers qubits that aren't I/O to a gate

    It takes two arguments and uses a dynamic programming approach to fill out a 2D array.

    Args:
        self (Layer) : this layer configurations object
        unsorted_pauli_ops (List[PauliOperator]) : a list of all the possible operators at this layer of the circuit
        
    Returns:
        void : fills self.pos_to_fill and self.carry_over_qubits
    """

    def check_qubits(self, unsorted_pauli_ops:List[PauliOperator]):
        self.pos_to_fill = defaultdict(list) # Hash map, where each key is a PauliOperator object 
        # and the key's associated value is the list of all the gate positions with non-identity I/O 
        # between the PauliOperator key and the layer to which we propagate.

        self.carry_over_qubits = [] # Represents the non-gate qubits of each PauliOperator object

        for i in range(len(unsorted_pauli_ops)): # For each valid configuration of our layer
            self.carry_over_qubits.append(copy.deepcopy(unsorted_pauli_ops[i].operator))
            for ind1, ind2 in self.gate_pos: # For each gate between this layer and its neighboring layer
                if (unsorted_pauli_ops[i].operator[ind1] != 'I' or unsorted_pauli_ops[i].operator[ind2] != 'I'): # non-identity output
                    self.pos_to_fill[unsorted_pauli_ops[i]].append((ind1,ind2)) # Adds gate to the list of positions 
                    # for layer i that require non-identity input
                    self.carry_over_qubits[i][ind1] = 'I'
                    self.carry_over_qubits[i][ind2] = 'I' # So when we compare the qubits that don't change,
                    # we don't also compare the ones that can change in a gate

    def group_sibs(self,unsorted_pauli_ops:List[PauliOperator]):
        sorted_pauli_ops = defaultdict(list) # Hash map of PauliOperators,
        # where each key (tuple of an entry from pos_to_fill and carry_over_qubits) maps to a list of 
        # all PauliOperators with a particular set of non-gate qubits and non-identity I/O gate positions

        for i in range(len(unsorted_pauli_ops)): # For each PauliOperator
            identifier = (tuple(self.pos_to_fill[unsorted_pauli_ops[i]]), tuple(self.carry_over_qubits[i]))
            
            sorted_pauli_ops[identifier].append(unsorted_pauli_ops[i]) # DefaultFict handles new keys
            # by initializing a new list for that key before trying to append
        
        if (self.backward):
            self.backward_rnp_sibs = sorted_pauli_ops
        else:
            self.forward_rnp_sibs = sorted_pauli_ops