#!pip install import-ipynb
import unittest
from typing import List, Tuple, DefaultDict
from collections import defaultdict
from pauli_operator import PauliOperator
from pauli_op_layer import PauliOpLayer
from pauli_path_trav import PauliPathTrav
from circuit_sim import CircuitSim
from prob_dist import ProbDist
from Brute_Force_RCS.circuit_utils import random_circuit, extract_gates_info

class TestProbDist(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.numQubits = 4
        self.depth = 2
        circuit = CircuitSim(self.numQubits, self.depth, [[(0,1), (2,3)], [(1,2)]]) # 1D
        self.bruteForceQC = random_circuit(self.numQubits, self.depth) # Qiskit Representation of a random circuit.
        gates = extract_gates_info(self.bruteForceQC)
        self.prob_dist = ProbDist(circuit, gates)

    def test_initialization(self):
        print()
        

if __name__ == '__main__':
    unittest.main()