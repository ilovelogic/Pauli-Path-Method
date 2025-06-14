# Fourier Coefficient and Marginal Sampling

This module implements the core calculation for the Fourier-based quasi-probability approximation \(\bar{q}(C, x)\) used to simulate noisy quantum circuits. It evaluates path-summed Fourier coefficients using legal Pauli operator sequences (Pauli paths) and supports marginal sampling using Lemma 9 from Aharonov et al. The module is  integrated with the Pauli path generator and supports both full output and partial (marginal) probability estimation.

---

## Table of Contents

- [Introduction](#introduction)
- [Architecture Overview](#architecture-overview)
- [Function Documentation](#function-documentation)
  - [compute_noisy_fourier](#compute_noisy_fourier)
  - [compute_marginal_noisy_fourier](#compute_marginal_noisy_fourier)
- [Supporting Functions](#supporting-functions)
- [Usage Example](#usage-example)
- [Equations](#equations)
- [Future Improvements](#future-improvements)

---

## Introduction

This module provides two main capabilities:

- Compute the quasi-probability \(\bar{q}(C, x)\) for a specific output string \(x\) using truncated Pauli paths.
- Compute marginal probabilities over a subset of qubits by reordering summation and applying Lemma 9 from Aharonov et al.

These features enable efficient simulation of noisy random quantum circuits, even under the presence of depolarizing noise.

---

## Architecture Overview

### Fourier Coefficient Computation

- Traverses the Pauli path trees (generated using sibling operators or XYZGenerations).
- At each layer, computes the trace transition amplitude \(\text{Tr}(s_i U_i s_{i-1} U_i^\dagger)\).
- Applies noise attenuation factor \((1 - \gamma)^{|s_i|}\) at each step.
- Final contribution is weighted by input and output overlaps (or partial overlap if computing marginals).

### Marginal Sampling

- Uses bit-by-bit sampling strategy:
  - At each step, compute conditional marginals \(\bar{q}(x_i = 0 \mid x_{<i})\) and \(\bar{q}(x_i = 1 \mid x_{<i})\)
  - Normalize, sample, and append to fixed bits.
- Based on the summation reordering enabled by Lemma 9.

---

## Function Documentation

### `compute_noisy_fourier(C, sib_op_heads, x, n, gamma)`

**Purpose:**  
Compute the Fourier coefficient sum \(\bar{q}(C, x)\) for a full output bitstring \(x\).

**Parameters:**

| Name          | Type         | Description |
|---------------|--------------|-------------|
| `C`           | list         | Preprocessed circuit layers with gates |
| `sib_op_heads`| list         | Pauli path tree roots |
| `x`           | str          | Output bitstring (e.g., "0110") |
| `n`           | int          | Number of qubits |
| `gamma`       | float        | Depolarizing noise rate |

**Returns:** `float` — Approximated output probability

---

### `compute_marginal_noisy_fourier(C, sib_op_heads, fixed_bits, n, gamma)`

**Purpose:**  
Compute the marginal probability \(\sum_{x: x_T = \text{fixed\_bits}} \bar{q}(C, x)\) where `fixed_bits` specifies the known qubit outcomes.

**Parameters:**

| Name          | Type         | Description |
|---------------|--------------|-------------|
| `C`           | list         | Preprocessed circuit layers |
| `sib_op_heads`| list         | Pauli path tree roots |
| `fixed_bits`  | dict         | Partial bitstring (e.g., `{0: '1', 2: '0'}`) |
| `n`           | int          | Number of qubits |
| `gamma`       | float        | Depolarizing noise rate |

**Returns:** `float` — Marginal output probability

---

## Supporting Functions

These helper functions support the main computations:

- `calculate_gate_transition_amplitude(sd, sd_minus_1, gate, qubits)`
- `calculate_layer_transition_amplitude(sd, sd_minus_1, layer_gates, n)`
- `calculate_input_overlap(s0)`
- `calculate_output_overlap(x, sd)`
- `calculate_partial_overlap(fixed_bits, sd)`
- `preprocess_circuit_gates(raw_gate_data, n)`

Each ensures proper normalization of Pauli matrices and compatibility with Qiskit’s endianness.

---

## Usage Example

```python
from Lemma_8.Pauli_Amplitude.fourier_coeff import compute_noisy_fourier, compute_marginal_noisy_fourier
from Lemma_8.path_generation import generate_pauli_path_tree
from Brute_Force_RCS.circuit_utils import extract_gates_info, random_circuit

# Generate circuit and extract gates
qc = random_circuit(num_qubits=4, depth=3)
gates = extract_gates_info(qc)
layers = preprocess_circuit_gates(gates, 4)

# Generate legal Pauli paths with truncation ℓ
sib_op_heads = generate_pauli_path_tree(num_qubits=4, depth=3, hamming_weight_cutoff=2)

# Compute full output probability
prob_x = compute_noisy_fourier(layers, sib_op_heads, x="0110", n=4, gamma=0.01)

# Compute marginal probability for qubit 0 = 1 and qubit 2 = 0
marginal_prob = compute_marginal_noisy_fourier(layers, sib_op_heads, fixed_bits={0: '1', 2: '0'}, n=4, gamma=0.01)
