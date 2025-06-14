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
## How Fourier Coefficient Calculation Works

Let’s walk through a concrete example showing how the Fourier coefficient is computed using Pauli path integrals and trace transition amplitudes.

### Problem Setup

- **Number of qubits:** 4  
- **Circuit depth:** 3  
- **Output bitstring (`x`)**: Last input to the function  
- **Pauli path sequence (`s`)**: `["IZIZ", "YIIX", "ZZIX", "ZZIZ"]`

**Circuit Description**:
Layer 1: CNOT(0,1), CZ(2,3)
Layer 2: SWAP(0,1)
Layer 3: CNOT(2,3)

markdown
Copy
Edit

We are computing the expression:  
`f(C, s, x) = Tr(|x⟩⟨x| ⋅ s_d) ⋅ ∏_{i=1}^d Tr(s_i ⋅ U_i ⋅ s_{i−1} ⋅ U_i†) ⋅ Tr(s₀ ⋅ |0ⁿ⟩⟨0ⁿ|)`

The total trace is a product of:
1. Input overlap
2. Layer-wise transition amplitudes
3. Output overlap (or marginal overlap)

### Breakdown of Layer Transition (Depth d = 3)

Suppose we are computing the transition amplitude for layer 3:  
Gate applied: `CNOT(2,3)`  
Current Pauli: `s_d = ZZIZ`  
Previous Pauli: `s_{d−1} = ZZIX`

**Step-by-step:**

1. **Extract sub-Paulis**:  
   - Gate acts on qubits (2,3)
   - So we extract characters at indices 2 and 3 from both `s_d` and `s_{d-1}`
     - `sd_sub = 'IZ'`
     - `sdm1_sub = 'IX'`

2. **Convert to normalized matrices**:  
   - These substrings are converted to `2x2` matrices using Qiskit’s `Pauli().to_matrix()` and normalized by dividing by `sqrt(2^len(substring))`

3. **Apply unitaries**:  
   - Apply the gate: `CNOT @ s_{d−1} @ CNOT†`
   - Then take the trace: `Tr(s_d ⋅ result)`

4. **Handle qubits not acted on by gates**:  
   - Indices (0,1) are idle/non-gate qubits, meaning that no gate is acting on them
   - So we also compute:  
     `Tr(sd_sub ⋅ sdm1_sub)` for the idle qubits:  
     `Tr(ZZ ⋅ ZZ)` which gives `2` (identity-like behavior), multiplied into the layer result

5. **Combine**:  
   - `Tr_total = Tr(IZ ⋅ CNOT ⋅ IX ⋅ CNOT†) * Tr(ZZ ⋅ I ⋅ ZZ ⋅ I)`

> Note: Anne’s path generator ensures we only get legal Pauli paths compatible with the circuit topology, so we already skip paths that would give zero contributions (e.g., those with `X/Y` on idle qubits).

### Relevant Functions That Handle This

| Function | Description |
|---------|-------------|
| `calculate_layer_transition_amplitude(sd, sd_minus_1, layer_gates, n_qubits)` | Orchestrates layer-wise transition amplitude computation |
| `calculate_gate_transition_amplitude(sd, sd_minus_1, gate, qubit_indices)` | Computes `Tr(s ⋅ U ⋅ s' ⋅ U†)` for active 2-qubit gates |
| `calculate_non_gate_transition_amplitude(sd, sd_minus_1, qubit_indices)` | Handles idle qubit transitions via `Tr(s ⋅ s')` |
| `extract_qubit_pauli(pauli_string, qubits)` | Utility to extract substrings from Pauli strings based on gate indices |

Each layer’s amplitude is multiplied recursively during Pauli path traversal (see `traverse_tree_with_noise`), and the final result is aggregated into the Fourier coefficient.

---

This modular decomposition helps keep the Fourier simulation fast and interpretable. Trace properties like `Tr(A ⊗


## Compute Noisy Fourier: `compute_noisy_fourier`

This function calculates the full Fourier-based output probability for a given bitstring `x` by recursively summing over all legal Pauli paths and computing:

f(C, s, x) = Tr(|x⟩⟨x| ⋅ s_d) ⋅ ∏ Tr(s_i ⋅ U_i ⋅ s_{i−1} ⋅ U_i†) ⋅ Tr(s₀ ⋅ |0ⁿ⟩⟨0ⁿ|)

### Inputs:

- `C`: List of circuit layers
- `sib_op_heads`: Pauli path trees (from Anne's generator)
- `x`: Bitstring string (e.g., `'0110'`)
- `n`: Number of qubits
- `gamma`: Depolarizing noise rate

- Traverses all valid Pauli paths using `traverse_tree_with_noise()`
- Applies depolarizing noise using factor `(1 - gamma)^|s|`
- Accumulates and sums all non-zero contributions

**Function used**:
```python
compute_noisy_fourier(C, path_roots, x, n, gamma)

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