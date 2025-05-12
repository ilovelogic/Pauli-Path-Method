import numpy as np
from qiskit.quantum_info import Pauli, Statevector
from qiskit import QuantumCircuit
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt

def reverse_qubit_indices(pauli_path):
    """Reverse qubit indices in a Pauli path to match Qiskit's little-endian convention."""
    return [op[::-1] for op in pauli_path]

def reverse_output_state(x):
    """Reverse output state string to match Qiskit's qubit ordering."""
    return x[::-1]


def normalize_pauli(pauli_string):
    """
    Normalize a Pauli string to match the normalized Pauli basis.
    Each single-qubit Pauli operator is divided by sqrt(2).
    """
    n = len(pauli_string)  # Number of qubits
    normalization_factor = 1 / np.sqrt(2**n)
    pauli_matrix = Pauli(pauli_string).to_matrix()
    
    return normalization_factor * pauli_matrix 


def extract_qubit_pauli(pauli_string, qubits):
    """
    Extracts the relevant sub-Pauli operator for specific qubits.
    
    Parameters:
        pauli_string (str): Full Pauli string (e.g., "IZIZ").
        qubits (list): Indices of qubits to extract (e.g., [0, 1]).
    
    Returns:
        str: Sub-Pauli operator for specified qubits.
    """
    # print(f"extract_qubit_pauli: pauli_string='{pauli_string}', qubits={qubits}")
    for q in qubits:
        if q >= len(pauli_string):
            print(f" Index {q} out of range for string of length {len(pauli_string)}")
    return ''.join([pauli_string[q] for q in qubits])

def calculate_gate_transition_amplitude(sd, sd_minus_1, gate, qubit_indices):
    """
    Calculate Tr(s_d U_gate s_{d-1} U_gate†) for a single 2-qubit gate.
    
    Parameters:
        sd (str): Current full Pauli operator (e.g., "IZIZ").
        sd_minus_1 (str): Previous full Pauli operator (e.g., "XXIX").
        gate (np.ndarray): Unitary matrix for a single 2-qubit gate.
        qubit_indices (list): Indices of qubits that the gate acts on.
    
    Returns:
        float: Transition amplitude for this gate.
    """
    # Extract relevant sub-Paulis
    sd_sub = extract_qubit_pauli(sd, qubit_indices)
    sd_minus_1_sub = extract_qubit_pauli(sd_minus_1, qubit_indices)
    
    # Convert to normalized matrices
    sd_matrix = Pauli(sd_sub).to_matrix() / np.sqrt(2**len(sd_sub))
    sd_minus_1_matrix = Pauli(sd_minus_1_sub).to_matrix() / np.sqrt(2**len(sd_minus_1_sub))
    
    # Compute transition amplitude
    gate_unchanged = gate
    gate_conj = gate.conj().transpose()
    transformed_sd_minus_1 = gate_unchanged @ sd_minus_1_matrix @ gate_conj
    return np.trace(sd_matrix @ transformed_sd_minus_1) #.real

def calculate_non_gate_transition_amplitude(sd, sd_minus_1, qubit_indices):
    """
    Calculate transition amplitude for qubits not acted upon by any gate.
    
    Parameters:
        sd (str): Current full Pauli operator.
        sd_minus_1 (str): Previous full Pauli operator.
        qubit_indices (list): Indices of qubits not acted upon by gates.
    
    Returns:
        float: Transition amplitude for non-gate qubits.
    """
    if not qubit_indices:  # If empty, return 1.0
        return 1.0
        
    # Extract relevant sub-Paulis
    sd_sub = extract_qubit_pauli(sd, qubit_indices)

    sd_minus_1_sub = extract_qubit_pauli(sd_minus_1, qubit_indices)
    
    # Convert to normalized matrices
    sd_matrix = Pauli(sd_sub).to_matrix() / np.sqrt(2**len(sd_sub))
    sd_minus_1_matrix = Pauli(sd_minus_1_sub).to_matrix() / np.sqrt(2**len(sd_minus_1_sub))
    
    # For non-gate qubits, the transition is just Tr(sd_sub · sd_minus_1_sub)
    return np.trace(sd_matrix @ sd_minus_1_matrix) #.real

def calculate_layer_transition_amplitude(sd, sd_minus_1, layer_gates, n_qubits):
    """
    Calculate transition amplitude for an entire layer of gates.
    
    Parameters:
        sd (str): Current full Pauli operator (e.g., "IZIZ").
        sd_minus_1 (str): Previous full Pauli operator (e.g., "XXIX").
        layer_gates (list of tuples): List of gates in the layer. Each tuple contains:
            - The unitary matrix for a 2-qubit gate.
            - The indices of qubits that the gate acts on.
        n_qubits (int): Total number of qubits in the system.
    
    Returns:
        float: Total transition amplitude for this layer.
    """
    layer_amplitude = 1.0
    
    # Track which qubits are acted upon by gates
    acted_qubits = set()
    for _, qubit_indices in layer_gates:
        acted_qubits.update(qubit_indices)
    
    # Calculate transition amplitude for each gate
    for gate, qubit_indices in layer_gates:
        gate_amplitude = calculate_gate_transition_amplitude(sd, sd_minus_1, gate, qubit_indices)
        layer_amplitude *= gate_amplitude
    
    # Calculate transition amplitude for qubits not acted upon by any gate
    non_gate_qubits = [i for i in range(n_qubits) if i not in acted_qubits]
    if non_gate_qubits:
        non_gate_amplitude = calculate_non_gate_transition_amplitude(sd, sd_minus_1, non_gate_qubits)
        layer_amplitude *= non_gate_amplitude
    
    return layer_amplitude

def calculate_input_overlap(s0):
    """
    Calculate Tr(s0 |0^n><0^n|).
    
    Parameters:
        s0 (str): Initial Pauli operator (e.g., "IZIZ").
    
    Returns:
        float: Input overlap.
    """
    # Check if s0 contains only I and Z (legal condition)
    if not all(op in ['I', 'Z'] for op in s0):
        return 0.0
        
    n = len(s0)  # Number of qubits
    
    # For a legal initial Pauli operator with only I and Z:
    # Count Z operators to determine sign
    z_count = s0.count('Z')
    
    # Normalization factor
    norm_factor = 1.0 / np.sqrt(2**n)
    
    # Return normalized overlap with sign based on Z count
    return norm_factor #* ((-1)**z_count)

def calculate_output_overlap(x, sd):
    """
    Calculate Tr(|x><x| s_d).
    
    Parameters:
        x (str): Output state as a binary string (e.g., "0000").
        sd (str): Final Pauli operator (e.g., "ZZII").
    
    Returns:
        float: Output overlap.
    """
    # Check if sd contains only I and Z (legal condition)
    if not all(op in ['I', 'Z'] for op in sd):
        return 0.0
    
    n = len(sd)  # Number of qubits
    
    # Normalization factor
    norm_factor = 1.0 / np.sqrt(2**n)
    
    # Calculate sign based on Z operators and corresponding bits in x
    sign = 1
    for i, op in enumerate(sd):
        if op == 'Z' and x[i] == '1':
            sign *= -1
            
    return norm_factor * sign


def compute_fourier_coefficient(C, s, x):
    """
    Compute f(C, s, x) for a given circuit C and Pauli path s.
    
    Parameters:
        C (list of list of tuples): Circuit as a list of layers,
                                    where each layer contains tuples of gates and their acting qubits.
                                    Example: [[(CNOT01, [0, 1]), (CNOT23, [2, 3])], ...]
        s (list of str): Pauli path as a list of strings (e.g., ["IZIZ", "YIXI", "XXIX", "ZZII"]).
        x (str): Output state as a binary string (e.g., "0000").
    
    Returns:
        float: Fourier coefficient f(C, s, x).
    """

    # Reverse qubit indices for Qiskit compatibility
    s_reversed = reverse_qubit_indices(s)
    x_reversed = reverse_output_state(x)
    n = len(s[0])  # Number of qubits
    d = len(s)-1     # Depth of the circuit
    
    # Check if path is legal (s0 and sd contain only I and Z)
    if not all(op in ['I', 'Z'] for op in s[0]) or not all(op in ['I', 'Z'] for op in s[-1]):
        return 0.0
    
    # Input overlap
    input_overlap = calculate_input_overlap(s_reversed[0])
    #input_overlap = calculate_input_overlap(s[0])
    if input_overlap == 0:
        return 0.0
    
    # Transition amplitudes
    transition_amplitude = 1.0
    for i in range(d):
        #layer_amplitude = calculate_layer_transition_amplitude(s[i+1], s[i], C[i], n)
        layer_amplitude = calculate_layer_transition_amplitude(
            s_reversed[i+1], 
            s_reversed[i],
            C[i],
            n
        )
        transition_amplitude *= layer_amplitude
        if transition_amplitude == 0:
            return 0.0
    
    # Output overlap
    output_overlap = calculate_output_overlap(x_reversed, s_reversed[-1])
    #output_overlap = calculate_output_overlap(x, s[-1])
    
    return input_overlap * transition_amplitude * output_overlap

# Preprocessing functions for taking in anne and jesus input 
def preprocess_circuit_gates(raw_gate_data, n):
    from collections import defaultdict
    # Reverse map: logical qubit i → tensor axis (bit) num_qubits - 1 - i
    #reversed_qubit_map = {q: n - 1 - q for q in range(n)
    layers = defaultdict(list)
    #for each 2qubit gate in a given layer
    for gate_matrix, qubits, layer in raw_gate_data:
        #mapped_qubits = [reversed_qubit_map[q] for q in qubits]

        reversed_qubits = [n - 1 - q for q in qubits]  # Reverse gate indices
        #print(f"Original qubits: {qubits}, Reversed: {reversed_qubits}")
        if len(qubits) == 2:
            #q0, q1 = reversed_qubits
            sorted_qubits = sorted(reversed_qubits)
            # For CNOT gates specifically
            if np.allclose(gate_matrix, np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])):
                # This transforms the matrix using basis state effects, not tensor reshaping
                new_matrix = np.zeros((4,4))
                # Preserve |00⟩ → |00⟩ and |10⟩ → |10⟩
                new_matrix[0,0] = 1
                new_matrix[2,2] = 1
                # Swap |01⟩ → |11⟩ and |11⟩ → |01⟩
                new_matrix[1,3] = 1
                new_matrix[3,1] = 1
                gate_matrix = new_matrix
            else:
                # For non-CNOT gates like RXX
                gate_matrix = np.reshape(gate_matrix, (2, 2, 2, 2))
                if reversed_qubits != sorted_qubits:
                    gate_matrix = gate_matrix.transpose((1, 0, 3, 2))
                gate_matrix = gate_matrix.reshape(4, 4)
        layers[layer].append((gate_matrix, reversed_qubits))
    return [layers[i] for i in sorted(layers)]

def preprocess_pauli_path(raw_path):
    return [''.join(layer) for layer in raw_path]

def compute_fourier_from_raw_inputs(raw_gate_data, raw_pauli_path, output_state, n):
    circuit_layers = preprocess_circuit_gates(raw_gate_data, n)  # Pass `n` here
    pauli_path_str = preprocess_pauli_path(raw_pauli_path)
    reversed_pauli_path = reverse_qubit_indices(pauli_path_str)
    reversed_output = reverse_output_state(output_state)
    #return compute_fourier_coefficient(circuit_layers, pauli_path_str, output_state)
    return compute_fourier_coefficient(circuit_layers, reversed_pauli_path, reversed_output)
