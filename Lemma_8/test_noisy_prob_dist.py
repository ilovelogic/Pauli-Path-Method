#!pip install import-ipynb
import unittest
from typing import List, Tuple, DefaultDict
from collections import defaultdict
from circuit_sim import CircuitSim
from Brute_Force_RCS import circuit_utils
from prob_dist import ProbDist
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


class TestNoisyProbDist(unittest.TestCase):
    #python -m Lemma_8.test_prob_dist
    @classmethod
    def setUpClass(self):
        
        self.numQubits = 4 # must be at least 3
        self.depth = 2
        
        self.C = circuit_utils.random_circuit(self.numQubits, self.depth)
        print()
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        print(self.C)

        self.bruteForceQC = self.C # Qiskit Representation of a random circuit.
        gates = circuit_utils.extract_gates_info(self.bruteForceQC)

        gate_pos = []

        #endian order indexing
        for i in range(len(gates)):
            layer_num = gates[i][2]
            if layer_num+1 > len(gate_pos):
                gate_pos.append([])
            a, b = gates[i][1]
            gate_pos[layer_num].append((self.numQubits - a - 1, self.numQubits - b - 1))
            #gate_pos[layer_num].append(gates[i][1])

        noise_rate = 0.001

        print("\n Noise rate: %s", noise_rate)
        print("Truncation parameter: (depth+1)*numQubits-2*1000000")

        start = time.time()

        circuit = CircuitSim(self.numQubits, (self.depth+1)*self.numQubits, gate_pos) # 1D, keeps all paths
        self.prob_dist = ProbDist(circuit, gates, self.numQubits, self.depth, self.bruteForceQC, noise_rate)

        end = time.time()
        duration = end - start

        print("\n \n Time taken for pauli prob dist generation: ")
        print(duration)

        print("prob dist: ")
        print(self.prob_dist.probs)
        return

    def test_stat_measures(self):
        return
        #self.assertEqual(0,self.prob_dist.tvd)
        #self.assertEqual(1,self.prob_dist.xeb)


if __name__ == '__main__':
    unittest.main()