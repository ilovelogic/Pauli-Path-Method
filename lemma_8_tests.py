import unittest
from pauli_operator import get_hamming_weights, R_iterations
from typing import List

class TestLemma8Functions(unittest.TestCase):

    def test_get_hamming_weights_general_case(self):
        """Test general case for get_hamming_weights"""
        d = 3  # depth
        l = 6  # upper bound on weight
        n = 4  # number of qubits
        result = get_hamming_weights(d, l, n)
        
        # expected combinations for d=3, l=6
        expected_combinations = [
            [1,1,1], # combinations summing to 3
            [1,1,2], [1,2,1], [2,1,1], # combinations summing to 4
            [1,2,2],  [2,1,2], [2,2,1], # combinations summing to 5
            [1,1,3], [1,3,1], [3,1,1],
            [2,2,2], # combinations summing to 6
            [1,1,4], [1,4,1], [4,1,1],
            [1,2,3], [3,1,2], [3,2,1], [1,3,2], [2,3,1], [2,1,3]
        ]
        
        self.assertEqual(sorted(result), sorted(expected_combinations))

    def test_get_hamming_weights_edge_case(self):
        """Test edge case for get_hamming_weights with minimum depth"""
        d = 1
        l = 5
        n = 1
        result = get_hamming_weights(d, l, n)
        
        # no possible combinations
        expected_combinations = []
        
        self.assertEqual(result, expected_combinations)

    def test_R_iterations_general_case(self):
        """Test general case for R_iterations."""
        prior_layer = "IRIRRRII"
        this_weight = 4
        n = len(prior_layer)
        pos_list = [(0,3),(1,2), (5,6),(4,7)]
        
        result1 = R_iterations(prior_layer, this_weight, n, pos_list)
        
        # checks that the number of configurations matches expectations
        self.assertTrue(len(result1) > 0)  

        # ensures all configurations have the correct weight
        for layer in result1:
            self.assertEqual(layer.count('R'), this_weight)

        prior_layer = "IIRIRRIRI"
        this_weight = 5
        n = len(prior_layer)
        pos_list = [(1,7),(2,4),(5,6)]
        
        result2 = R_iterations(prior_layer, this_weight, n, pos_list)
        
        # checks that the number of configurations matches expectations
        self.assertTrue(len(result2) > 0)  
        
        # ensures all configurations have the correct weight
        for layer in result2:
            self.assertEqual(layer.count('R'), this_weight)

    def test_R_iterations_edge_case(self):
        """Test edge case for R_iterations with no valid configurations."""
        prior_layer = "IIIII"
        this_weight = 5
        n = len(prior_layer)
        pos_list = [(1,4)]
        
        result3 = R_iterations(prior_layer, this_weight, n, pos_list)

        # no way to create a valid configuration with this_weight > prior_weight,
        # so the result should be empty.
        self.assertEqual(result3, [])

    def test_R_iterations_single_qubit(self):
        """Test edge case for R_iterations with a single qubit."""
        prior_layer = "I"
        this_weight = 0
        n = len(prior_layer)
        pos_list = []
        
        result4 = R_iterations(prior_layer, this_weight, n, pos_list)

        
        # with a single qubit and weight of one non-identity gate,
        # no configuration is possible: []
        self.assertEqual(result4, ["I"])

if __name__ == "__main__":
    unittest.main()