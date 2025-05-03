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

import pdb

from qiskit import QuantumCircuit 

#python -m pdb your_script.py
#break exception IndexError

class TestProbDist(unittest.TestCase):
    #python -m Lemma_8.test_prob_dist
    @classmethod
    def setUpClass(self):

        return
        self.numQubits = 3 # must be at least 3
        self.depth = 1

        self.C = circuit_utils.random_circuit(self.numQubits, self.depth)

        self.bruteForceQC = self.C # Qiskit Representation of a random circuit.
        gates = circuit_utils.extract_gates_info(self.bruteForceQC)

        gate_pos = []
        for i in range(len(gates)):
            layer_num = gates[i][2]
            if layer_num+1 > len(gate_pos):
                gate_pos.append([])
            gate_pos[layer_num].append(gates[i][1])
        print(gate_pos)

        circuit = CircuitSim(self.numQubits, (self.depth+1)*self.numQubits, gate_pos) # 1D, keeps all paths
        
        self.prob_dist = ProbDist(circuit, gates, self.bruteForceQC)

    #def test_stat_measures(self):
        # self.assertEqual(0,self.prob_dist.tvd)
        self.assertEqual(1,self.prob_dist.xeb)
        #self.assertEqual(0,self.prob_dist.tvd)
        print(self.prob_dist.tvd)
        self.assertEqual(1,self.prob_dist.xeb)
        print(self.prob_dist.xeb)

    def test_no_depth(self):
        self.n_qubits_no_depth(3)

    def n_qubits_no_depth(self, n:int):
        self.numQubits = n # must be at least 3
        self.depth = 1 # we will only have identities, so it's effectively a depthless circuit

        self.QC = QuantumCircuit(n)
        for i in range(n):
            self.QC.id(i)

        self.bruteForceQC = self.QC
        gates = circuit_utils.extract_gates_info(self.QC)
        noise_rate = 0

        s_list = [[list(p), list(p)] for p in product('IZ', repeat=n)]

        total_prob = 0
        self.probs = DefaultDict(float) # hash function mapping outcomes to their probabilities
        for i in range(1 << self.numQubits):
            x = format(i, f'0{self.numQubits}b') # possible outcome of the circuit, represented as a string of 1's and 0's
            self.probs[x] = 0

            for s in s_list:
                ham_weight = 1 # total number of non-identity Paulis in s
                # each non-identity Pauli is affected by the depolarizing noise
                # E(ρ) := (1 − γ)ρ + γ(I/2)Tr(ρ)
                #pdb.set_trace()
                self.probs[x] += ((1-noise_rate)**ham_weight)*compute_fourier_from_raw_inputs(gates, s, x)
            print(f'p({x}) = {self.probs[x]}')
            total_prob += self.probs[x]
        print(f'Total probability sum = {total_prob}') # sum should be 1

if __name__ == '__main__':
    unittest.main()