# pip install "qiskit-aer>=0.11.0"
from collections import defaultdict
from typing import List, DefaultDict
from pauli_operator import PauliOperator
from circuit_sim import CircuitSim
import numpy as np
from Brute_Force_RCS.evaluation_utils import total_variation_distance, calculate_true_distribution, compute_xeb
from Brute_Force_RCS.circuit_utils import  complete_distribution, generate_emp_distribution
from Pauli_Amplitude.pauli_amplitude import compute_fourier_coefficient

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
        self.C = gates # list of tuples, containing the layer of each gate, the matrix, and the qubit indicices its acting on
        self.probs = DefaultDict(float) # made hash function
        n = circuit.num_qubits
        for i in range(1 << n):
            x= format(i, f'0{n}b') # possible outcome of the circuit, represented as a string of 1's and 0's
            self.probs[x] = 0
            for s in self.s_list:
                self.probs[x] += compute_fourier_coefficient(self.C, s, x)
        self.calc_TVD()
        self.calc_linearXEB()

    # ------------------------------------------------------------------------------
    # TVD of pauli prob dist and true dist

    def calc_TVD(self):
      #TVD of true distribution and pauli probability distribution
    
      trueDist = calculate_true_distribution(self.bruteForceQC)
      # trueDist assumes that we can access the qiskit representation of whatever 1D
      # circuit we generated.
      # Im assuming the self class can contain the 1d circuit

      full_prob_dist = complete_distribution(self.prob_dist)
      # full prob dist just ensures that every possible basis state is present in the
      # probability outcome to work with my TVD function.
      ideal_TVD = total_variation_distance(trueDist, full_prob_dist)
      return ideal_TVD

    def calc_linearXEB(self):
      #XEB of true distribution and pauli probability distribution

      trueDist = calculate_true_distribution(self.bruteForceQC)
      full_prob_dist = complete_distribution(self.prob_dist)
      xeb = compute_xeb(trueDist, full_prob_dist)
      return xeb

    def pauli_ops_to_strs(self, xyz_pauli_paths:List[List[List[str]]]):
        self.s_list = [[] for _ in range(len(xyz_pauli_paths))]
        for i in range(len(xyz_pauli_paths)):
            for pauli_op in xyz_pauli_paths[i]:
                self.s_list[i].append(pauli_op.operator)
