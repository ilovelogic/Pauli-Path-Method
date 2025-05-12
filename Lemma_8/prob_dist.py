# pip install "qiskit-aer>=0.11.0"
from collections import defaultdict
from typing import List, DefaultDict
from pauli_operator import PauliOperator
from circuit_sim import CircuitSim
import numpy as np
import pdb
from Brute_Force_RCS.evaluation_utils import total_variation_distance, calculate_true_distribution, compute_xeb
from Brute_Force_RCS.circuit_utils import  complete_distribution, generate_emp_distribution
from Pauli_Amplitude.pauli_amplitude import compute_fourier_from_raw_inputs
from qiskit import circuit
import itertools
import copy


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
        self.num_qubs = num_qubs
        self.depth = depth

        #test on this one, right now the values aren't looking right 
        self.pauli_ops_to_strs(circuit_sim.xyz_pauli_paths) # initializes self.s_list, which contains all pauli paths

        self.C = gates # list of tuples, containing the layer of each gate, the matrix, and the qubit indicices its acting on
        self.probs = DefaultDict(float) # hash function mapping outcomes to their probabilities
        self.n = circuit_sim.num_qubits
        #not going to the right probability states for this one 
        self.bruteForceQC = QC

        self.s_list = self.brute_force_paths()
        
        self.calc_noisy_prob_dist(noise_rate)
    
        self.calc_TVD()
        self.calc_linearXEB()

        


    # Algorithm 1 from the rcs paper
    def calc_noisy_prob_dist(self, noise_rate:float):
      total_prob = 0
      for i in range(1 << self.n):
        x = format(i, f'0{self.n}b') # possible outcome of the circuit, represented as a string of 1's and 0's
        
        self.probs[x] = 0

        for s in self.s_list:
          #print(s)
          ham_weight = self.get_hamming_weight(s) # total number of non-identity Paulis in s
          # each non-identity Pauli is affected by the depolarizing noise
          # E(ρ) := (1 − γ)ρ + γ(I/2)Tr(ρ)
          fourier_coeff = compute_fourier_from_raw_inputs(self.C, s, x, self.n)
          self.probs[x] += ((1-noise_rate)**ham_weight)*fourier_coeff
          #print out the pauli paths 
          #if (abs(fourier_coeff) > 1/(10**10)):
             #print(f'{x} and {s}, amplitude = {fourier_coeff}')
        
        #printing the output state from erika's code and total prob
        print(f'p({x}) = {self.probs[x]}')
        total_prob += self.probs[x]
      print(f'Total probability sum = {total_prob}') # sum should be 1

    # Determines number of non-identity Paulis in a given legal Pauli path
    def get_hamming_weight(self,path:List[List[str]]):
      hamming_weight = 0
      for pauli_op in path:
         for pauli in pauli_op:
            if pauli != 'I':
               hamming_weight += 1
      return hamming_weight
    
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
      inner_list_length = self.num_qubs
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
          #print(f"\n=== Layer {i} Pauli Paths ===")
          for pauli_op in xyz_pauli_paths[i]:
              #pauli_str = ''.join(pauli_op.operator)
              self.s_list[i].append(pauli_op.operator)
              #print(f"Path {j}: {pauli_str}")


        # accounting for the fact that we excluded the all I's case from our path generation
        self.s_list[len(xyz_pauli_paths)] = [["I" for _ in range(len(xyz_pauli_paths[0][0].operator))] for _ in range(len(xyz_pauli_paths[0]))]    
        #for lil_list in self.s_list:
           #print(lil_list)

        for i in range(len(self.s_list)):
            path = copy.deepcopy(self.s_list[i])
            self.s_list[i] = [''.join(layer) for layer in path]
