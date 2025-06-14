from pauli_amplitude import compute_fourier_coefficient, compute_fourier_from_raw_inputs
import unittest
from Lemma_8.circuit_sim import CircuitSim
from Lemma_8.get_prob_dist import GetProbDist
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
        gate_01 = np.array([
            [ 0.52276082+0.0361418j ,  0.152548  +0.19747164j, 0.33408134+0.67642091j,  0.24221037+0.18795612j],
            [-0.00945917+0.18095772j, -0.17053029-0.7138475j , -0.1078335 +0.25962698j,  0.42519478-0.41070813j],
            [-0.2382395 +0.53491495j,  0.40765612+0.37185391j, 0.17755533+0.00582369j,  0.05121   -0.56432864j],
            [-0.52144587+0.27889399j, -0.04551292-0.30421203j, 0.51916506+0.22520659j, -0.28771874+0.39072178j]
        ])

        gate_23 = np.array([
            [ 0.16933099+0.49518196j, -0.17648941+0.09440534j, -0.25311353-0.72416522j, -0.10429706-0.29445095j],
            [ 0.0954568 -0.3989199j ,  0.13016058-0.06222619j, -0.25997905-0.16294595j,  0.76582755-0.36097745j],
            [ 0.51149066+0.43312816j,  0.4806417 -0.2285176j ,  0.28233598+0.28847464j,  0.06998835-0.31577044j],
            [-0.32031373-0.07777638j,  0.71827872-0.37411789j, -0.20561242-0.33491806j, -0.19753231+0.2049289j ]
        ])

        raw_gate_data = [
            (gate_01, (0, 1), 0),
            (gate_23, (2, 3), 0),
            (gate_01, (0, 1), 1),
            (gate_23, (2, 3), 2),
        ]

        raw_pauli_path = [
            ['Z', 'Z', 'I', 'Z'],
            ['I', 'X', 'Y', 'I'],
            ['I', 'I', 'Z', 'I'],
            ['I', 'I', 'I', 'Z']
        ]

        output_state = "0000"

        result1 = compute_fourier_from_raw_inputs(raw_gate_data, raw_pauli_path, output_state)
        print(f" Fourier coefficient with raw inputs f(C, s, x) = {result1}")
        self.assertEqual(result1, 0.0)
    #Fourier coefficient f(C, s, x) = 0.0 result should be 0 

    def test_compute_fourier_coefficient2(self):
        CZ = np.array([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, -1]
        ])

        # Raw gate data across 3 layers
        raw_gate_data_CZ = [
            (CZ, (0, 1), 0),
            (CZ, (2, 3), 0),
            (CZ, (0, 1), 1),
            (CZ, (2, 3), 1),
            (CZ, (0, 1), 2),
            (CZ, (2, 3), 2),
        ]

        # Raw Pauli path as list of lists
        raw_pauli_path_Z = [
            ['Z', 'Z', 'Z', 'Z'],
            ['Z', 'Z', 'Z', 'Z'],
            ['Z', 'Z', 'Z', 'Z'],
            ['Z', 'Z', 'Z', 'Z']
        ]

        # Output state
        output_state = "0000"

        # Run the computation
        result2 = compute_fourier_from_raw_inputs(raw_gate_data_CZ, raw_pauli_path_Z, output_state)
        print(f" Fourier coefficient with raw inputs f(C, s, x) = {result2}")
        self.assertEqual(result2, 0.0625)





if __name__ == '__main__':
    unittest.main()