import numpy as np
from qiskit.quantum_info import Pauli, Statevector
from qiskit import QuantumCircuit
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt
import copy
from typing import Set

import pdb # my debugging bestie


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
    #z_count = s0.count('Z')
    
    # Normalization factor
    norm_factor = 1.0 / np.sqrt(2**n)
    
    # Return normalized overlap with sign based on Z count
    return norm_factor #* ((-1)**z_count)

def calculate_partial_overlap(fixed_bits, sd):
    """
    Tr(sd ⋅ (⨂_{i ∈ T} |x_i⟩⟨x_i| ⊗ ⨂_{j ∉ T} I)) from Lemma 9.
    Returns the marginal measurement term for fixed qubits.
    """
    n = len(sd)
    k = len(fixed_bits)

    #corrected_bits = {n - 1 - i: b for i, b in fixed_bits.items()}

    #making sure it's a legal s_d path, only Is and Zs
    if not all(p in ['I', 'Z'] for p in sd):
        return 0.0
    

    sign = 1
    for i, p in enumerate(sd):
        if i in fixed_bits:
            if p == 'Z' and fixed_bits[i] == '1':
                sign *= -1
        elif i not in fixed_bits and p != 'I':
            return 0.0


    return sign * (1 / np.sqrt(2 ** n)) * (2 ** (n - k))



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

def compute_marginal_noisy_fourier(C, sib_op_heads, fixed_bits, n, gamma):
    """
    Computes ∑_{x ∈ {0,1}^n: x_T = fixed_bits} q̄(C, x) as per Lemma 9.
    """

    fourier_coeffs_for_paths = []
    #visited_roots = set()  

    for root in sib_op_heads:
        '''
        if id(root) in visited_roots:
            continue
        visited_roots.add(id(root)) 
        '''
        traverse_tree_with_noise(
            root,
            fourier_coeffs_for_paths,
            1.0,
            [],
            -1,
            C,
            x=None,
            n=n,
            gamma=gamma,
            fixed_bits=fixed_bits
        )
    #print(f"[DEBUG] Total paths contributing to marginal {fixed_bits}: {len(fourier_coeffs_for_paths)}")

    return sum(fourier_coeffs_for_paths)


def compute_noisy_fourier(C, sib_op_heads, x, n, gamma):
    #pdb.set_trace()
    """
    Compute the total Fourier coefficient f(C, s, x) by traversing the SiblingOps tree.
    
    Parameters:
        C (list): Preprocessed circuit as layers of (unitary, [qubits]) tuples.
        sib_op_heads (List[SiblingOps]): Root nodes of Pauli path trees.
        x (str): Output bitstring (e.g., "0000").
        n (int): Number of qubits.
        
    Returns:
        float: Sum of Fourier coefficients over all legal Pauli paths.
    """

    #x = reverse_output_state(x) # reverse the output state to match Qiskit
    
    total = 0.0 # probability for this state
    fourier_coeffs_for_paths = [] # all fourier coefficients for paths ending at this state
    #pdb.set_trace()

    for root in sib_op_heads:
        #build_list(root, fourier_coeffs_for_paths, 1, 
                                    #[], -1, C, x, n, gamma)
        traverse_tree_with_noise(root, fourier_coeffs_for_paths, 1, [], -1, C, x, n, gamma)

    for fourier_coeff in fourier_coeffs_for_paths:
        #print(fourier_coeff)
        total += fourier_coeff
    #print(f'fourier_coeffs_for_paths = {len(fourier_coeffs_for_paths)}')
    return total

def build_list(sib, fourier_coeffs_for_paths, cur_fourier, op_list, index, C, x, n, gamma):
    """
    Recursively traverse a SiblingOps tree to accumulate Fourier coefficient contributions.
    
    Parameters:
        sib_op (SiblingOps): Current node.
        path_so_far (List[str]): List of Pauli strings (s0 to sd-1 so far).
        C (list): Circuit layers.
        x (str): Output bitstring.
        total (List[float]): Single-element list used to accumulate total.
        n (int): Number of qubits.
    """
    for op in sib.pauli_ops:

        cur_op = op.operator #[::-1] # reverses the operator to match Qiskit

        this_op_list = copy.deepcopy(op_list)
                
        this_op_list.append(cur_op)
  

        #ham_weight = sum(p != 'I' for p in cur_op) # accounting for noise
        #branched_cur_fourier *= (1 - gamma) ** ham_weight
        
        if sib.next_sibs is None:
            print(f'index = {index}, op_list = {this_op_list}')
            fourier_coeffs_for_paths.append(this_op_list)
        else: 
            for next_sib in sib.next_sibs:
                build_list(next_sib, fourier_coeffs_for_paths, cur_fourier, this_op_list, index+1, C, x, n, gamma)

def traverse_tree_with_noise(sib, fourier_coeffs_for_paths, cur_fourier, prev_op, index, C, x, n, gamma, 
                             fixed_bits=None):
    """
    Recursively traverse a SiblingOps tree to accumulate Fourier coefficient contributions.
    
    Parameters:
        sib_op (SiblingOps): Current node.
        path_so_far (List[str]): List of Pauli strings (s0 to sd-1 so far).
        C (list): Circuit layers.
        x (str): Output bitstring.
        total (List[float]): Single-element list used to accumulate total.
        n (int): Number of qubits.
    """

    for op in sib.pauli_ops:


        cur_op = op.operator #[::-1] # reverses the operator to match Qiskit

        if prev_op == []:
            layer_amplitude = calculate_input_overlap(cur_op)

        else:
            layer_amplitude = calculate_layer_transition_amplitude(
                cur_op, 
                prev_op,
                C[index],
                n
            )
        
        #print(f"[DEBUG] amplitude {layer_amplitude} at depth {index+1} with ops: {prev_op} → {cur_op}")
 
        # ask how the reversal comes into play
        branched_cur_fourier = cur_fourier * layer_amplitude

        # each non-identity Pauli is affected by the depolarizing noise
        # E(ρ) := (1 − γ)ρ + γ(I/2)Tr(ρ)
        ham_weight = sum(p != 'I' for p in cur_op) # accounting for noise
        branched_cur_fourier *= (1 - gamma) ** ham_weight
        
        cur_op_tuple = tuple(cur_op)
        key = (index, cur_op_tuple)
        #if key in visited:
        #    return
        #visited.add(key)
        #print(f"[DEBUG] depth {index+1} | cur_op = {cur_op}, prev_op = {prev_op}")
        #MAX_DEPTH = len(C)

        if sib.next_sibs is None:
            #print(f"[DEBUG] FINAL layer | cur_op = {cur_op}")

            if fixed_bits is not None:
                final_fourier = branched_cur_fourier * calculate_partial_overlap(fixed_bits, cur_op)
            else:
                final_fourier = branched_cur_fourier * calculate_output_overlap(x, cur_op)

            if final_fourier != 0:
                fourier_coeffs_for_paths.append(final_fourier)
            #continue

        else:
            for next_sib in sib.next_sibs:
                #print(f"[RECURSE] Going deeper: index = {index+1}, cur_op = {cur_op}")
                traverse_tree_with_noise(
                        next_sib,
                        fourier_coeffs_for_paths,
                        branched_cur_fourier,
                        cur_op,
                        index + 1,
                        C,
                        x,
                        n,
                        gamma,
                        fixed_bits,
                    )


        ''' 
        if sib.next_sibs is None:
            final_fourier = branched_cur_fourier * calculate_output_overlap(x, cur_op)
            if final_fourier != 0:
                fourier_coeffs_for_paths.append(final_fourier)
                #print(f'outcome = {x}, prev = {prev_op}, cur = {cur_op}, amplitude = {final_fourier}')
            #elif final_fourier == 0:
                #print(f'0 = {x}, prev = {prev_op}, cur = {cur_op}, index = {index}')
        else: 
            for next_sib in sib.next_sibs:
                traverse_tree_with_noise(next_sib, fourier_coeffs_for_paths, branched_cur_fourier, cur_op, index+1, C, x, n, gamma)
     '''