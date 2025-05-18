# pip install "qiskit-aer>=0.11.0"
from collections import defaultdict
from typing import List, DefaultDict
from pauli_operator import PauliOperator
from circuit_sim import CircuitSim
import numpy as np
import pdb
from Brute_Force_RCS.evaluation_utils import total_variation_distance, calculate_true_distribution, compute_xeb
from Brute_Force_RCS.circuit_utils import  complete_distribution, generate_emp_distribution
from Pauli_Amplitude.og_pauli_amp import compute_fourier_from_raw_inputs, preprocess_circuit_gates
from Pauli_Amplitude.edited_pauli_amp import compute_noisy_fourier
from qiskit import circuit
import itertools
# pip3 install memory-profiler requests
from memory_profiler import profile
import requests


class ProbDist:
    """
    
    """
    def __init__(self, circuit_sim:CircuitSim,gates:List, num_qubs:int, depth:int, QC:circuit, noise_rate:float=0):

        '''
        circuit (CircuitSim): A fully initiated CircuitSim object based on our circuit architecture
        gates (List[(int,np.ndarray,tuple)): the first item in the tuple is the layer it is in (0 indexing),
        the second item is the gate matrix, and the third item is a tuple of the gate positions
        QC (circuit): quantum circuit generated using Qiskit
        noise_rate: single-qubit depolarizing noise (γ in the research paper)
        '''
        self.depth = depth
        self.n = num_qubs

        #test on this one, right now the values aren't looking right 
        self.pauli_ops_to_strs(circuit_sim.xyz_pauli_paths) # initializes self.s_list, which contains all pauli paths

        #self.C = gates
        self.C = preprocess_circuit_gates(gates, self.n) # list of tuples, containing the layer of each gate, the matrix, and the qubit indicices its acting on
        self.probs = DefaultDict(float) # hash function mapping outcomes to their probabilities
        

        # tree roots for Pauli path traversal
        self.sib_op_heads = circuit_sim.sib_op_heads

        self.bruteForceQC = QC
        
        self.calc_noisy_prob_dist(noise_rate)
    
        self.calc_TVD()
        self.calc_linearXEB()

        


    # Algorithm 1 from the rcs paper
    @profile
    def calc_noisy_prob_dist(self, noise_rate:float):

      self.other_probs = DefaultDict(float) # hash function mapping outcomes to their probabilities
      for i in range(1 << self.n):
        x = format(i, f'0{self.n}b') # possible outcome of the circuit, represented as a string of 1's and 0's
   
        self.probs[x] = compute_noisy_fourier(self.C, self.sib_op_heads, x, self.n, noise_rate)
        self.probs[x] += compute_fourier_from_raw_inputs(self.C, 
                        [["I" for _ in range(self.n)] for _ in range(len(self.C)+1)], x, self.n)
        print(f'p({x}) = {self.probs[x]}')

    # ------------------------------------------------------------------------------
    # TVD of pauli prob dist and true dist

    def calc_TVD(self):
      #TVD of true distribution and pauli probability distribution

      trueDist = calculate_true_distribution(self.bruteForceQC)
      # trueDist assumes that we can access the qiskit representation of whatever 1D
      # circuit we generated.
      # Im assuming the self class can contain the 1d circuit

      full_prob_dist = complete_distribution(self.probs,self.n)
      # full prob dist just ensures that every possible basis state is present in the
      # probability outcome to work with my TVD function.
      self.tvd = total_variation_distance(trueDist, full_prob_dist) # replace with outs


    def calc_linearXEB(self):
      #XEB of true distribution and pauli probability distribution

      trueDist = calculate_true_distribution(self.bruteForceQC)
      print(trueDist)
      full_prob_dist = complete_distribution(self.probs,self.n)
      self.xeb = compute_xeb(trueDist, full_prob_dist, self.n)


    def brute_force_paths(self):

      elements = ['X', 'Y', 'Z', 'I']
      inner_list_length = self.n
      outer_list_length = self.depth+1

      # generates all possible inner lists of length 3 (for 3 qubits)
      all_inner_lists = list(itertools.product(elements, repeat=inner_list_length))

      # builds all possible lists of 3 inner lists (for the 3 Pauli operators in a path)
      all_outer_lists = list(itertools.product(all_inner_lists, repeat=outer_list_length))

      # converts tuples to lists since that is the format Erika's code uses
      all_outer_lists_as_lists = [ [list(inner) for inner in outer] for outer in all_outer_lists]

      return all_outer_lists_as_lists


    def pauli_ops_to_strs(self, xyz_pauli_paths:List[List[List[str]]]):
        
        self.s_list = [[] for _ in range(len(xyz_pauli_paths)+1)]
        for i in range(len(xyz_pauli_paths)):
          for pauli_op in xyz_pauli_paths[i]:
              self.s_list[i].append(pauli_op.operator)


        # accounting for the fact that we excluded the all I's case from our path generation
        self.s_list[len(xyz_pauli_paths)] = [["I" for _ in range(len(xyz_pauli_paths[0][0].operator))] for _ in range(len(xyz_pauli_paths[0]))]

        #for lil_list in self.s_list:
           #print(lil_list)
