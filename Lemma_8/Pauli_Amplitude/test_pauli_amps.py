from pauli_amplitude import compute_fourier_coefficient
import unittest
from Lemma_8.circuit_sim import CircuitSim
from Lemma_8.prob_dist import ProbDist
from Brute_Force_RCS.circuit_utils import random_circuit, extract_gates_info
from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer
from qiskit import assemble
from qiskit.visualization import plot_histogram
from qiskit.quantum_info import Pauli
import numpy as np
import matplotlib.pyplot as plt

class TestProbDist(unittest.TestCase):
    def test_compute_fourier_coefficient(self):
        # Define gates
        CNOT = np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 1, 0]])
        
        CZ = np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, -1]])
        
        SWAP = np.array([[1, 0, 0, 0],
                        [0, 0, 1, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 1]])
        
        # Define circuit
        circuit = [
            [(CNOT, [0, 1]), (CZ, [2, 3])],    # Layer 1
            [(SWAP, [0, 1])],                  # Layer 2
            [(CNOT, [2, 3])]                   # Layer 3
        ]
        
        # Define Pauli path
        pauli_path = ["IZIZ", "YIIX", "ZZIX", "ZZIZ"]

        
        # Define output state
        output_state = "0000"
        
        # Calculate Fourier coefficient
        result = compute_fourier_coefficient(circuit, pauli_path, output_state)
        print(f"Fourier coefficient f(C, s, x) = {result}")
    #Fourier coefficient f(C, s, x) = 0.0
        

if __name__ == '__main__':
    unittest.main()