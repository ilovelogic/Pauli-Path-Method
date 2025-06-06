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
        
        truncation_param = 2
        l = (self.depth+1)*self.numQubits - truncation_param

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

        print()
        print(f'Noise rate: {noise_rate}')
        print(f'Truncation parameter: {truncation_param}')

        start = time.time()

        circuit = CircuitSim(self.numQubits, l, gate_pos) # 1D, keeps all paths
        self.prob_dist = ProbDist(circuit, gates, self.numQubits, self.depth, self.bruteForceQC, noise_rate)

        end = time.time()
        duration = end - start


        print(f'\n \nTime taken for Pauli probability distribution generation:')
        print(duration)
       
        return

    def test_stat_measures(self):
        #print(f'XEB={self.prob_dist.xeb}')
        #print(f'TVD={self.prob_dist.tvd}')
        return
        #self.assertEqual(0,self.prob_dist.tvd)
        #self.assertEqual(1,self.prob_dist.xeb)


if __name__ == '__main__':
    unittest.main()