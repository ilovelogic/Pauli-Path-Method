# Noisy Output Probability Distribution Calculation

The `NoisyProbDist` class uses the `GetProbDist` class to derive noisy quantum circuit output distributions using the Pauli path method, as described in [Aharonov et al.](https://arxiv.org/pdf/2211.03999). This builds on the code described in the README files of `Lemma 8`, `Fourier Coefficient`, and `Brute Force RCS`. It combines them to provide the next step in the simulation pipeline: computing output probabilities and benchmarking them against brute-force simulation.

---

## Table of Contents

- [Introduction](#introduction)
- [Architecture Overview](#architecture-overview)
- [Class Documentation](#class-documentation)
  - [GetProbDist](#getprobdist)
  - [NoisyProbDist](#noisyprobdist)
- [Usage Example](#usage-example)
- [Output Metrics](#output-metrics)
- [References](#references)

---

## Introduction

The `NoisyProbDist` and `GetProbDist` classes implement the core of the Pauli path method for simulating noisy random circuit sampling (RCS). `NoisyProbDist` takes in general circuit specifications and uses them to build a random quantum circuit. It then instantiates a `CircuitSim` to obtain all possible legal Pauli paths in line with the upperbound on Hamming weight. 

Then it feeds this information into `GetProbDist`, which calculates the output probability distribution, incorporating depolarizing noise, and compares the results to brute-force Qiskit simulations using metrics such as TVD (Total Variation Distance), XEB (Linear Cross-Entropy Benchmark), and fidelity.

---

## Architecture Overview

1. **NoisyProbDist**:  
   - Generates a random circuit, extracts gate data, builds the Pauli path structure, and runs `GetProbDist` with user-specified parameters.
   - Reports computation time and quality metrics.

2. **GetProbDist**:  
    - Consumes a precomputed set of Pauli paths and circuit/gate information, which are sent in as arguments to the constructor.
    - Computes the output probability distribution for all bitstrings using the Pauli path integral, with optional depolarizing noise.
    - Calculates TVD, XEB, and fidelity against brute-force Qiskit simulation.

---

## Class Documentation

### ProbDist

**Purpose:**  
Computes the noisy or noiseless output probability distribution of a quantum circuit using the Pauli path method, and compares the results to brute-force simulation.

**Initialization:**
ProbDist(
circuit_sim: CircuitSim,
gates: List,
num_qubs: int,
depth: int,
QC: QuantumCircuit,
noise_rate: float = 0
)

**Parameters:**
| Parameter    | Type           | Description                                                |
|--------------|----------------|------------------------------------------------------------|
| circuit_sim  | CircuitSim     | Precomputed legal Pauli paths for the circuit              |
| gates        | list           | List of gate tuples (layer, matrix, qubit indices)         |
| num_qubs     | int            | Number of qubits                                           |
| depth        | int            | Circuit depth                                              |
| QC           | QuantumCircuit | Qiskit circuit object for brute-force simulation           |
| noise_rate   | float          | Depolarizing noise parameter (γ), default 0 (noiseless)    |

**Key Methods:**
- `calc_noisy_prob_dist()`:  
  Computes the output probability for every bitstring using the Pauli path integral, including noise if specified.
- `calc_TVD()`:  
  Calculates the Total Variation Distance between the computed distribution and the brute-force Qiskit distribution (noisy or noiseless as appropriate).
- `calc_linearXEB()`:  
  Computes the Linear Cross-Entropy Benchmark between the computed and true distributions.
- `calc_fidelity()`:  
  Computes the classical fidelity between the computed and true distributions.

**Attributes:**
- `probs`: `DefaultDict[str, float]` — Output probability for each bitstring.
- `tvd`: `float` — Total Variation Distance.
- `xeb`: `float` — Linear XEB.
- `fidelity`: `float` — Classical fidelity.

---

### NoisyProbDist

**Purpose:**  
Generates a random quantum circuit, extracts gates and Pauli paths, and computes the noisy output probability distribution using `ProbDist`.

**Initialization:**
NoisyProbDist(
num_qubits: int,
depth: int,
truncation_param: int,
noise_rate: float = 0.001
)


**Parameters:**
| Parameter         | Type | Description                                        |
|-------------------|------|----------------------------------------------------|
| num_qubits        | int  | Number of qubits (≥3)                              |
| depth             | int  | Circuit depth                                      |
| truncation_param  | int  | Hamming weight cutoff for Pauli paths              |
| noise_rate        | float| Per-qubit depolarizing noise (default: 0.001)      |

**Attributes:**
- `n`: `int` — The predetermined number of qubits for our circuit.
- `d`: `int` — The specified circut depth.
- `l`: `int` — The upperbound on Hamming weight.
- `bruteForceQC`: `QuantumCircuit` — Qiskit Representation of a random circuit with number of qubits `n` and gate depth `d`.
- `prob_dist`: `GetProbDist` — The probability distribution object for the circuit.
- `duration`: `float` — Time taken (seconds) for probability distribution generation.
- `bruteForceQC`: `QuantumCircuit` — The randomly generated Qiskit circuit.

---

>## Usage Example
    >from Brute_Force_RCS.circuit_utils import random_circuit, extract_gates_info\
    >from Lemma_8.get_prob_dist import GetProbDist, NoisyProbDist\
    >\
    >dist = NoisyProbDist(4, 3, 2, 0.001) \# 4-qubit, depth-3, truncation=2, noise=0.001\
    >\
    > \#Access probability distribution and metrics\
    >print(f"TVD: {dist.prob_dist.tvd:.4f}")\
    >print(f"XEB: {dist.prob_dist.xeb:.4f}")\
    >print(f"Fidelity: {dist.prob_dist.fidelity:.4f}")\
    >print(f"Time taken: {dist.duration:.2f} seconds")



---

## Output Metrics

- **Probability Distribution:**  
  Output probability for every computational basis state.
- **Total Variation Distance (TVD):**  
  Distance between computed and brute-force distributions.
- **Linear XEB:**  
  Cross-entropy between computed and true distributions.
- **Fidelity:**  
  Classical fidelity between distributions.
- **Computation Time:**  
  Wall-clock time for probability distribution generation.

---

## References

- [Aharonov et al., "Polynomial-Time Classical Simulation of Noisy Random Circuit Sampling"](https://arxiv.org/pdf/2211.03999)