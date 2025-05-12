#!pip install import-ipynb
import unittest
from typing import List, Tuple, DefaultDict
from collections import defaultdict
from circuit_sim import CircuitSim
from Brute_Force_RCS import circuit_utils
from prob_dist import ProbDist
from Pauli_Amplitude.pauli_amplitude import compute_fourier_from_raw_inputs
from qiskit import circuit
from itertools import product
from Brute_Force_RCS.evaluation_utils import total_variation_distance, calculate_true_distribution, compute_xeb
from Brute_Force_RCS.circuit_utils import  complete_distribution, generate_emp_distribution
import math
import numpy as np

import pdb

from qiskit import QuantumCircuit 


class TestProbDist(unittest.TestCase):
    #python -m Lemma_8.test_prob_dist
    @classmethod
    def setUpClass(self):
        
        self.numQubits = 3 # must be at least 3
        self.depth = 1

        #self.C = QuantumCircuit(self.numQubits)
        # making two-qubit  H⊗I matrix
        #H = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]])  # Hadamard
        #I = np.eye(2)                                     # Identity
        #H_tensor_I = np.kron(H, I)  # Applies H to qubit 0, I to qubit 1

        # applies it to qubits 0 and 1
        #self.C.unitary(H_tensor_I, [2, 1], label='H⊗I')
        #self.C.unitary(H_tensor_I, [1, 0], label='H⊗I')
 
        #self.C.cx(1,2)
        #self.C.rxx(math.pi / 2,1,2)
    
        #print(self.C)
        #for i in range(1,self.numQubits-1,1):
            #self.C.cx(i,i+1)
        
        self.C = circuit_utils.random_circuit(self.numQubits, self.depth)

        self.bruteForceQC = self.C # Qiskit Representation of a random circuit.
        gates = circuit_utils.extract_gates_info(self.bruteForceQC)

        gate_pos = []

        #endian order indexing
        for i in range(len(gates)):
            layer_num = gates[i][2]
            if layer_num+1 > len(gate_pos):
                gate_pos.append([])
            reversed_qubits = [self.numQubits - 1 - q for q in gates[i][1]]
            gate_pos[layer_num].append(tuple(reversed_qubits))
            #a, b = gates[i][1]
            #gate_pos[layer_num].append((self.numQubits - a - 1, self.numQubits - b - 1))
            #gate_pos[layer_num].append(gates[i][1])

        #circuit = CircuitSim(self.numQubits, (self.depth+1)*self.numQubits, gate_pos) # 1D, keeps all paths
        circuit = CircuitSim(self.numQubits, self.depth+1, gate_pos) # 1D, keeps all paths
        print("YAY")
        
        self.prob_dist = ProbDist(circuit, gates, self.numQubits,self.depth, self.bruteForceQC)
        return

    def test_stat_measures(self):
        return
        self.assertEqual(1,self.prob_dist.xeb)
        self.assertEqual(0,self.prob_dist.tvd)



    def test_no_depth(self):
        return
        self.n_qubits_no_depth(3)

    def n_qubits_no_depth(self, n:int):
        return
        self.numQubits = n # must be at least 3
        self.depth = 1 # we will only have identities, so it's effectively a depthless circuit

        self.QC = QuantumCircuit(n)
        for i in range(n):
            self.QC.h(i)

        print(self.QC)

        self.bruteForceQC = self.QC
        gates = circuit_utils.extract_gates_info(self.QC)
        noise_rate = 0

        s_list = [[list(p), list(p)] for p in product('IZ', repeat=n)]

        total_prob = 0
        probs = DefaultDict(float) # hash function mapping outcomes to their probabilities
        for i in range(1 << self.numQubits):
            x = format(i, f'0{self.numQubits}b') # possible outcome of the circuit, represented as a string of 1's and 0's
            probs[x] = 0

            for s in s_list:
                ham_weight = 1 # total number of non-identity Paulis in s
                # each non-identity Pauli is affected by the depolarizing noise
                # E(ρ) := (1 − γ)ρ + γ(I/2)Tr(ρ)

                probs[x] += ((1-noise_rate)**ham_weight)*compute_fourier_from_raw_inputs(gates, s, x)
            print(f'p({x}) = {probs[x]}')
            total_prob += probs[x]
        print(f'Total probability sum = {total_prob}') # sum should be 1

        trueDist = calculate_true_distribution(self.bruteForceQC)
        print(trueDist)

        full_prob_dist = complete_distribution(probs,self.numQubits)

        tvd = total_variation_distance(trueDist, full_prob_dist) # replace with outs
        self.assertEqual(0,tvd)
        

if __name__ == '__main__':
    unittest.main()