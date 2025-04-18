from collections import defaultdict
from typing import List, DefaultDict
from pauli_operator import PauliOperator
from circuit_sim import CircuitSim
import numpy as np

class ProbDist:
    """
    
    """
    def __init__(self, circuit:CircuitSim,gates:List):
        '''
        circuit (CircuitSim): A fully initiated CircuitSim object based on our circuit architecture
        gates (List[(int,np.ndarray,tuple)): the first item in the tuple is the layer it is in (0 indexing),
        the second item is the gate matrix, and the third item is a tuple of the gate positions
        '''
        self.pauli_ops_to_strs(circuit.xyz_pauli_paths) # initializes self.s_list, which contains all pauli paths
        self.C = gates
        self.probs = DefaultDict(int) # made hash function
        n = circuit.num_qubits
        for i in range(1 << n):
            x= format(i, f'0{n}b') # possible outcome of the circuit, represented as a string of 1's and 0's
            self.probs[x] = 0
            for s in self.s_list:
                self.probs[x] += compute_fourier_coefficient(self.C, s, x)

        #for x in list[all out comes]:
            #self.probs.append(0)
            #for s in self.s_list:
                #self.probs[x] += compute_fourier_coefficient(self.C, s, x)

    def pauli_ops_to_strs(self, xyz_pauli_paths:List[List[List[str]]]):
        self.s_list = [[] for _ in range(len(xyz_pauli_paths))]
        for i in range(len(xyz_pauli_paths)):
            for pauli_op in xyz_pauli_paths[i]:
                self.s_list[i].append(pauli_op.operator)
