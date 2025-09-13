from typing import List, Tuple, DefaultDict
from collections import defaultdict
from Path_Generation.circuit_sim import CircuitSim
from Brute_Force_RCS import circuit_utils
from Prob_Calc.get_prob_dist import GetProbDist
from qiskit import circuit
from itertools import product
from Brute_Force_RCS.evaluation_utils import total_variation_distance, calculate_true_distribution, compute_xeb
from Brute_Force_RCS.circuit_utils import  complete_distribution, generate_emp_distribution
import math
import numpy as np
import warnings
import time

import pdb

from qiskit import QuantumCircuit 


class NoisyProbDist:
  """
  Generates a noisy probability distribution for a random quantum circuit.
  The circuit is represented as a QuantumCircuit object from Qiskit.
  """

  def __init__(self, num_qubits:int, depth:int, truncation_param:int, noise_rate=0.001):
    self.n = num_qubits # must be at least 3
    self.d = depth
    
    self.l = (self.d+1)*self.n - truncation_param

    self.bruteForceQC = circuit_utils.random_circuit(self.n, self.d) # Qiskit Representation of a random circuit

    print()
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    print(self.bruteForceQC)

    # Stores a list of tuples, where each tuple is (gate matrix, qubits applied to, layer number)
    # List is ordered so that gates on earlier layers appear first
    gates = circuit_utils.extract_gates_info(self.bruteForceQC)

    gate_pos = []

    for i in range(len(gates)):
        layer_num = gates[i][2]
        while layer_num+1 > len(gate_pos): # if gate_pos doesn't have a list for this layer yet
          gate_pos.append([])
        a, b = gates[i][1] # qubit positions where the gate is applied
        gate_pos[layer_num].append((self.n - a - 1, self.n - b - 1)) # reverse qubit order to be consistent with circuit_sim.py


    noise_rate = 0.001 # expected amount of noise per qubit per layer in real-world circuits

    print()
    print(f'Noise rate: {noise_rate}')
    print(f'Truncation parameter: {truncation_param}')

    start = time.time()

    circuit = CircuitSim(self.n, self.l, gate_pos) # 1D architecuture
    self.prob_dist = GetProbDist(circuit, gates, self.n, self.d, self.bruteForceQC, noise_rate)

    end = time.time()
    self.duration = end - start


    print(f'\n \nTime taken for Pauli probability distribution generation:')
    print(self.duration)