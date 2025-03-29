import unittest
from typing import List, Tuple, DefaultDict
from collections import defaultdict
from pauli_operator import PauliOperator
from pauli_op_layer import PauliOpLayer
from pauli_path_traversal import PauliPathTraversal
from circuit_sim import CircuitSim

class TestPauliPathTraversal(unittest.TestCase):
    def setUp(self):
        self.path = PauliPathTraversal(4, [1, 1, 1], [[(0, 1)], [(2, 3)]])

    def test_initialization(self):
        self.assertEqual(self.path.num_qubits, 3)
        self.assertEqual(self.path.depth, 2)
        self.assertEqual(self.path.weight_combo, [1, 1, 1])
        self.assertEqual(self.path.gate_pos, [[(0, 1)], [(1, 2)]])

    def test_build_min_configs(self):
        self.path.build_min_configs()
        self.assertIsInstance(self.path.layers, list)
        self.assertTrue(all(isinstance(layer, PauliOpLayer) for layer in self.path.layers))

    def test_unsorted_min_layer_ops(self):
        ops = self.path.unsorted_min_layer_ops(1)
        self.assertIsInstance(ops, list)
        self.assertTrue(all(isinstance(op, PauliOperator) for op in ops))

    def test_min_backward(self):
        min_layers = [PauliOpLayer([(0, 1)], 1)]
        self.path.min_backward(min_layers, 1)
        self.assertIsInstance(min_layers[0].backward_sibs, defaultdict)

    def test_min_forward(self):
        min_layers = [PauliOpLayer([(0, 1)], 0)]
        min_layer_ops = [PauliOperator(['I', 'R', 'I'])]
        self.path.min_forward(min_layers, min_layer_ops, 0)
        self.assertIsInstance(min_layers[0].forward_sibs, defaultdict)

    def test_propagate_next(self):
        all_sibs = defaultdict(list)
        all_sibs[(('R', 'I'), ('I',))] = [PauliOperator(['R', 'I', 'I'])]
        pos_to_fill = defaultdict(list)
        pos_to_fill[PauliOperator(['R', 'I', 'I'])] = [(0, 1)]
        next_ops = self.path.propagate_next(all_sibs, pos_to_fill, 0, 0)
        self.assertTrue(all(isinstance(op, PauliOperator) for op in next_ops.values))

class TestCircuit(unittest.TestCase):
    def setUp(self):
        self.circuit = CircuitSim(4, 10, [[(0, 1),(2,3)],[(1,2)]]) # brickwork 1D

    def test_initialization(self):
        self.assertEqual(self.circuit.num_qubits, 4)
        self.assertEqual(self.circuit.num_op_layers, 3)
        self.assertEqual(self.circuit.max_weight, 10)
        self.assertEqual(self.circuit.gate_pos, [[(0, 1),(2,3)],[(1,2)]])

    def test_init_pauli_paths(self):
        self.circuit.init_pauli_paths()
        self.assertIsInstance(self.circuit.pauli_paths, list)
        self.assertTrue(all(isinstance(path, PauliPathTraversal) for path in self.circuit.pauli_paths))

    def test_enumerate_weights(self):
        self.assertIsInstance(self.circuit.weight_combos, list)
        self.assertTrue(all(isinstance(combo, list) for combo in self.circuit.weight_combos))
        self.assertTrue(all(sum(combo) <= 10 for combo in self.circuit.weight_combos))
        self.assertEqual(len(self.circuit.weight_combos), len(set(tuple(combo) for combo in self.circuit.weight_combos)))

        # We want weight_combos to be the set of all lists where the ith element in each list represents 
        # the Hamming weight assigned to the ith Pauli operator in a legal Pauli path. 
        # Legal Pauli path Requirements:
        #   1) For each gate between Pauli operators, the I/O Paulis are either both 'II', or both not 'II'
        #   2) The first and last Pauli operators are only comprised of 'I's and 'Z's
        # Note that requirement 2 does not affect our weight lists. 
        # Requirement 1 means that the only way to decrease Hamming weight of one Pauli operator
        # respective to the prior Pauli operator is to have a gate have input 'RR' and output 'IR' or 'RI'.
        # Likewise, the sole way to increase a Pauli operator's Hamming weight respective to the preceeding 
        # operator's Hamming weight is to have a gate between these recieve 'IR' or 'RI' and output 'RR'. 
        # Non-gate Paulis are carried over to the next Pauli operator, so they do not increase or decrease the weight.
        true_weight_combos = {(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2), (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2), 
                            (1, 2, 3), (2, 2, 3), (2, 3, 2), (2, 3, 3), (3, 2, 1), (3, 2, 2), (3, 2, 3), (3, 3, 2), 
                            (3, 3, 3), (2, 3, 4), (2, 4, 3), (2, 4, 4),  (3, 3, 4), (3, 4, 3), (4, 2, 1), (4, 2, 2), 
                            (4, 2, 3), (4, 3, 2), (4, 3, 3)}
        self.assertEqual(true_weight_combos, set(tuple(combo) for combo in self.circuit.weight_combos))

if __name__ == '__main__':
    unittest.main()