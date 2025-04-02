import unittest
from typing import List, Tuple, DefaultDict
from collections import defaultdict
from pauli_operator import PauliOperator
from pauli_op_layer import PauliOpLayer
from pauli_path_trav import PauliPathTrav
from circuit_sim import CircuitSim

class TestCircuits(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.circuit = CircuitSim(4, 10, [[(0, 1),(2,3)],[(1,2)]]) # brickwork 1D

    def test_initialization(self):
        self.assertEqual(self.circuit.num_qubits, 4)
        self.assertEqual(self.circuit.num_op_layers, 3)
        self.assertEqual(self.circuit.max_weight, 10)
        self.assertEqual(self.circuit.gate_pos, [[(0, 1),(2,3)],[(1,2)]])

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

    def test_init_pauli_paths(self):
        self.circuit.init_pauli_paths()
        self.assertIsInstance(self.circuit.pauli_path_travs, list)
        self.assertTrue(all(isinstance(path, PauliPathTrav) for path in self.circuit.pauli_path_travs))

    def test_travs_to_list(self):
        pauli_path_set = set() # Set we will use to check for duplicate Pauli paths
        num_pauli_ops = 0 # Counter so we can compare the number of Pauli paths
        # stored in our list of lists of Pauli paths to the number of Pauli paths in the set
        for i in range (len(self.circuit.pauli_path_list)):
            this_weight_combo = self.circuit.weight_combos[i] # The weight combo assigned to every Pauli path in this list
            for pauli_path in self.circuit.pauli_path_list[i]:
                tuple_pauli = tuple(tuple(pauli_path_op.operator) for pauli_path_op in pauli_path) # Sets can store tuples, not lists
                pauli_path_set.add(tuple_pauli) # Only adds the Pauli path if it is not already in the set
                num_pauli_ops += 1
                # Checks if each Pauli operator in our Pauli path has the weight specified by its weight combo
                for i in range (len(pauli_path)):
                    pauli_op_weight = 0
                    for pauli in pauli_path[i].operator:
                        if pauli == 'R' or pauli == 'P':
                            pauli_op_weight += 1
                    if (pauli_op_weight != this_weight_combo[i]):
                        for pauli_op in pauli_path:
                            print(pauli_op.operator)
                        print(this_weight_combo)
                        print()
        self.assertEqual(len(pauli_path_set), num_pauli_ops) # If there are no duplicates, every element in 
        # the list of lists would have been added to the set, making these two values equal
    

if __name__ == '__main__':
    unittest.main()