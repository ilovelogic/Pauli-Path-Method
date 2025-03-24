from collections import defaultdict
import copy
from typing import DefaultDict, List
from pauli_operator import PauliOperator

class Layer:

    def __init__(self, gate_pos:List[tuple]=None, backward:int=-1,pauli_ops:DefaultDict[tuple, List[PauliOperator]]=None):
        """
        Initiate the forward and backward lists of operators
        for propagating forward and backward at a given depth of the circuit
        """
        
        self.backward = backward # if backward = 1, propagating backward
        # otherwise, propagating forward

        self.gate_pos = gate_pos # gate positions between this Layer and:
        # if (backward): the previous Layer
        # else: the next Layer

        if (pauli_ops is not None and self.gate_pos is not None and self.backward != -1):
            layer_ops = []

            # builds list of all PauliOperators at this depth of the circuit
            # turns DefaultDict into a 1D list
            for sibs in pauli_ops.values():
                for i in range(len(sibs)):
                    layer_ops.append(sibs[i])

            self.check_qubits(layer_ops)

            if (self.backward):
                self.forward_sibs = pauli_ops # the input list in this case comes sorted according to
                # the lists that the Layer at the next depth propagated to when it propagated backward
            else: 
                self.backward_sibs = pauli_ops # the input list in this case comes sorted according to
                # the lists that the Layer at the previous depth propagated to when it propagated forward
            
            self.find_sibs(layer_ops) # need to sort this list according to which ops at this depth propagate together 
            # either backward (if self.backward) or forward (else)

    """
    This function determines which gate positions of layers have non-identity I/O,
    and keeps track of each layers qubits that aren't I/O to a gate

    It takes two arguments and uses a dynamic programming approach to fill out a 2D array.

    Args:
        self (Layer) : this layer configurations object
        layers (List[PauliOperator]) : a list of all the possible layers at this depth of the circuit
        
    Returns:
        void : fills self.pos_to_fill and self.carry_over_qubits
    """

    def check_qubits(self, layers:List[PauliOperator]):
        self.pos_to_fill = defaultdict(list) # each layer's gate positions that require non-identity input
        # stored as a list of lists where first index is for layer and second is for gate

        self.carry_over_qubits = []

        for i in range(len(layers)): # for each valid configuration of our layer
            self.carry_over_qubits.append(copy.deepcopy(layers[i].operator))
            for ind1, ind2 in self.gate_pos: # for each gate between this layer and its neighboring layer
                if (layers[i].operator[ind1] == 'R' or layers[i].operator[ind2] == 'R'): # non-identity output
                    self.pos_to_fill[layers[i]].append((ind1,ind2)) # adds gate to the list of positions 
                    # for layer i that require non-identity input
                    self.carry_over_qubits[i][ind1] = 'I'
                    self.carry_over_qubits[i][ind2] = 'I' # so when we compare the qubits that don't change,
                    # we don't also compare the ones that can change in a gate


    def find_sibs(self,unsorted_pauli_ops:List[PauliOperator]):
        sorted_pauli_ops = defaultdict(list) # hash map of PauliOperators
        # where each key maps to a list of all PauliOperators 
        # with a particular set of non-gate qubits and non-identity input gate positions

        for i in range(len(unsorted_pauli_ops)): # for each PauliOperator
            identifier = (tuple(self.pos_to_fill[i]), tuple(self.carry_over_qubits[i]))
            # appends new items to the list associated with the identifier
            sorted_pauli_ops[identifier].append(unsorted_pauli_ops[i]) # defaultdict handles new keys
            # by initializing a new list for that key before trying to append43
        
        if (self.backward):
            self.backward_sibs = sorted_pauli_ops
        else:
            self.forward_sibs = sorted_pauli_ops