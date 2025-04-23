#!pip install import-ipynb
import unittest
from typing import List, Tuple, DefaultDict
from collections import defaultdict
from circuit_sim import CircuitSim
from Brute_Force_RCS import circuit_utils
from prob_dist import ProbDist

class TestProbDist(unittest.TestCase):
    #python -m Lemma_8.test_prob_dist
    @classmethod
    def setUpClass(self):
        self.numQubits = 4
        self.depth = 4
        
        self.C = circuit_utils.random_circuit(self.numQubits, self.depth)

        circuit = CircuitSim(self.numQubits, self.numQubits, [[(0,1), (2,3)], [(1,2)]]) # 1D, keeps all paths
        self.bruteForceQC = circuit_utils.random_circuit(self.numQubits, self.depth) # Qiskit Representation of a random circuit.
        gates = circuit_utils.extract_gates_info(self.bruteForceQC)
        print(gates)
        self.prob_dist = ProbDist(circuit, gates, self.bruteForceQC)

    def test_stat_measures(self):
        self.assertEqual(1,self.prob_dist.tvd)
        self.assertEqual(0,self.prob_dist.xeb)
        

if __name__ == '__main__':
    unittest.main()