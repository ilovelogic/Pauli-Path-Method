# Fourier Coefficient Calculation

The Fourier coefficient calculation module implements the core mathematical engine of the Pauli path method, which allows classically estimating noisy output probabilities from a quantum circuit. The implementation is based on the Fourier decomposition of quantum processes described in Aharonov et al., and it supports both full bitstring probability computations and marginal sampling over partial outputs.

---

## Table of Contents
- [Introduction](#introduction)
- [Key Equation](#key-equation)
- [Function Overview](#function-overview)
- [Usage Example](#usage-example)
- [Marginal Sampling](#marginal-sampling)
- [Performance Considerations](#performance-considerations)

---

## Introduction

This module enables the computation of noisy output probabilities `q̄(C, x)` for a given output bitstring `x`, or marginal probabilities over subsets of qubits, using a truncated Fourier expansion. The computation avoids summing over all `2^n` output strings and instead aggregates legal Pauli paths `s = (s_0, ..., s_d)` with bounded Hamming weight.

The implementation supports two use cases:
- `compute_noisy_fourier(C, heads, x, n, gamma)` — full bitstring probability
- `compute_marginal_noisy_fourier(C, heads, fixed_bits, n, gamma)` — marginal over a partial assignment

---

## Key Equation

The Fourier coefficient `f(C, s, x)` for a Pauli path `s` and bitstring `x` is computed as:

f(C, s, x) = Tr(|x⟩⟨x| ⋅ s_d) ⋅ ∏{i=1}^d Tr(s_i ⋅ U_i ⋅ s{i−1} ⋅ U_i†) ⋅ Tr(s₀ ⋅ |0ⁿ⟩⟨0ⁿ|)

To account for noise, each path is scaled by `(1 - γ)^|s|`, where `|s|` is the total number of non-identity Paulis.

---

## Function Overview

### `compute_noisy_fourier(C, heads, x, n, gamma)`
- **Purpose**: Computes `q̄(C, x)` for a specific bitstring `x`
- **Parameters**:
  - `C`: Circuit (as list of gate layers)
  - `heads`: Root nodes of the Pauli path trees
  - `x`: Output bitstring (e.g., `'0010'`)
  - `n`: Number of qubits
  - `gamma`: Depolarizing noise rate
- **Returns**: Approximate probability of observing bitstring `x`

### `compute_marginal_noisy_fourier(C, heads, fixed_bits, n, gamma)`
- **Purpose**: Computes marginal probability `∑_{x: x_T = fixed_bits} q̄(C, x)`
- **Parameters**:
  - `fixed_bits`: Dictionary `{qubit_index: '0' or '1'}` specifying a partial bitstring
- **Returns**: Marginal probability over all bitstrings consistent with `fixed_bits`

### Helper Functions
- `calculate_input_overlap(s0)`: Computes `Tr(s₀ ⋅ |0ⁿ⟩⟨0ⁿ|)`
- `calculate_output_overlap(x, s_d)`: Computes `Tr(|x⟩⟨x| ⋅ s_d)`
- `calculate_partial_overlap(fixed_bits, s_d)`: Computes marginal trace term per Lemma 9
- `calculate_layer_transition_amplitude(...)`: Handles product of traces across layers

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