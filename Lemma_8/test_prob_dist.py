#!pip install import-ipynb
import unittest
from typing import List, Tuple, DefaultDict
from collections import defaultdict
from circuit_sim import CircuitSim
from Brute_Force_RCS import circuit_utils
from prob_dist import ProbDist
import pdb
#python -m pdb your_script.py
#break exception IndexError

class TestProbDist(unittest.TestCase):
    #python -m Lemma_8.test_prob_dist
    @classmethod
    def setUpClass(self):
        self.numQubits = 4
        self.depth = 2
        
        self.C = circuit_utils.random_circuit(self.numQubits, self.depth)

        self.bruteForceQC = circuit_utils.random_circuit(self.numQubits, self.depth) # Qiskit Representation of a random circuit.
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

    def test_stat_measures(self):
        # self.assertEqual(0,self.prob_dist.tvd)
        self.assertEqual(1,self.prob_dist.xeb)
        #self.assertEqual(0,self.prob_dist.tvd)
        print(self.prob_dist.tvd)
        self.assertEqual(1,self.prob_dist.xeb)
        print(self.prob_dist.xeb)
        

if __name__ == '__main__':
    unittest.main()