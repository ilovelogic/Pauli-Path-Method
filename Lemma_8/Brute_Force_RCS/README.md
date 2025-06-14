# Documentation: Quantum Circuit Functions for Noisy RCS Simulation

This documentation explains each function in the provided code, emphasizing its role in simulating random quantum circuits as described in the paper by Aharonov et al. It highlights non-trivial behaviors such as measurement handling and their relevance to the simulation and analysis of noisy random circuit sampling (RCS).

---

## `random_circuit(num_qubits: int, depth=None)`

**Purpose:**  
Generates a Qiskit quantum circuit with a 1D brickwork architecture, where each two-qubit gate is a Haar-random unitary.

**Relation to Paper:**  
Constructs circuits that match the random circuit model discussed in the paper, essential for analyzing anti-concentration and noise effects. The use of Haar-random gates supports theoretical assumptions for anti-concentration and gate-set orthogonality.

**Notes/Kinks:**
- Defaults `depth` to ⌊log₂(n)⌋ if unspecified, matching the regime for expected anti-concentration.
- Measurement is **not** added by default; append manually if needed.

---

## `create_2d_brickwork_circuit(rows: int, cols: int, depth: int)`

**Purpose:**  
Builds a 2D brickwork quantum circuit over a specified grid with Haar-random two-qubit gates in a staggered pattern.

**Relation to Paper:**  
Implements a general architecture where anti-concentration can still hold, aligning with claims that results extend beyond 1D circuits.

---

## `extract_gates_info(qc)`

**Purpose:**  
Extracts matrix representations, qubit indices, and layer metadata for each gate.

**Relation to Paper:**  
Used in Pauli basis analysis and verifying circuit structure, which supports the paper’s path integral and Fourier framework.

---

## `run_ideal_simulation(qc: QuantumCircuit, shots=100)`

**Purpose:**  
Simulates the circuit on a noise-free backend and returns measurement outcome counts.

**Relation to Paper:**  
Provides the ideal RCS distribution `p_ideal(C, x)` for comparison with noisy outputs.

---

## `create_noise_model(depolarizing_param: float)`

**Purpose:**  
Creates a depolarizing noise model affecting all two-qubit gates.

**Relation to Paper:**  
Implements local depolarizing noise per gate, as specified in the noisy RCS model.

**Notes/Kinks:**
- Only valid for `0 < γ ≤ 1`.

---

## `run_noisy_simulation(circuit: QuantumCircuit, noise_model: NoiseModel, shots=100)`

**Purpose:**  
Simulates a quantum circuit under a noise model, returning measurement results.

**Relation to Paper:**  
Generates noisy distribution `p(C, x)`, which the classical algorithm is designed to approximate.

**Notes/Kinks:**
- Adds measurement before simulation and removes it afterward to support Qiskit's sampling requirements and enable repeated manipulation.

---

## `complete_distribution(partial_dist: dict, num_qubits: int) -> dict`

**Purpose:**  
Fills in missing computational basis states with zero-probability entries.

**Relation to Paper:**  
Ensures empirical distributions are defined over the full basis for total variation distance and statistical indistinguishability checks.

---

## `count_to_distribution(counts: dict[str, int], num_qubits, shots=100)`

**Purpose:**  
Normalizes raw counts into a full probability distribution over basis states.

**Relation to Paper:**  
Enables comparison of empirical vs. theoretical distributions, supporting evaluation of classical vs. quantum/noisy RCS outputs.

---

## `generate_emp_distribution(qc: QuantumCircuit, shots: int, noise=None, depth=None)`

**Purpose:**  
Runs the circuit and returns a normalized empirical distribution, optionally under noise.

**Relation to Paper:**  
Automates the generation of ideal or noisy RCS samples for benchmarking and statistical tests like total variation distance or XEB.

**Notes/Kinks:**
- Internally decides between ideal/noisy simulation based on the `noise` parameter.

---

## `reverse_keys(qiskit_basis_dist: dict[str, float]) -> dict[str, float]`

**Purpose:**  
Reverses bitstring keys to match endianness conventions.

**Relation to Paper:**  
Ensures Qiskit results align with theoretical conventions, which is crucial for correct statistical comparison.

---

## References

- Dorit Aharonov et al., *"A polynomial-time classical algorithm for noisy random circuit sampling"*, arXiv:2211.03999




# Evaluation Utilities for Noisy RCS Simulations

This documentation explains evaluation utilities for analyzing noisy Random Circuit Sampling (RCS) simulations. The metrics implemented (TVD, XEB, fidelity) are used to benchmark empirical data against theoretical predictions, as discussed in the paper:

**Reference:** Dorit Aharonov et al., *"A polynomial-time classical algorithm for noisy random circuit sampling"*, arXiv:2211.03999

---

## `calculate_true_distribution(qc)`

**Purpose:**  
Computes ideal (noiseless) output probabilities using statevector simulation.

**Relation to Paper:**  
Provides baseline distribution `p(C, x)` for distribution comparisons. 

**Notes/Kinks:**
- Must remove measurements before simulation (measurement treated as final step).
- Bitstrings are reversed to match theoretical conventions (Fig. 1).

---

## `total_variation_distance(dist1, dist2)`

**Purpose:**  
Computes  
![image](https://github.com/user-attachments/assets/efa5dd4a-18cb-4259-b47f-9823f053b435)


between two distributions.

**Relation to Paper:**  
Main metric for simulation accuracy (Theorem 1). Noisy RCS distributions are proven to be ϵ-close in TVD to classical outputs.

**Notes/Kinks:**
- Assumes complete distributions; use `complete_distribution` if necessary.
- Direct implementation of Eq. (26).

---

## `compute_xeb(dist1, dist2, num_qubits)`

**Purpose:**  
Calculates cross-entropy benchmarking (XEB) score:  


**Relation to Paper:**  
Implements XEB (linear Cross Entropy Benchmark)

![image](https://github.com/user-attachments/assets/19fe1a4a-9d67-43fe-981b-754ce4f5738a)


**Notes/Kinks:**
- Requires reversed keys for correct bitstring order.
- Matches the statistical test refuted by Theorem 1.

---

## `tvd_truedist_empdist(num_qubits, noise, shots)`

**Purpose:**  
Computes total variation distance between ideal and noisy empirical distributions.

**Relation to Paper:**  
Operationalizes Theorem 1’s result about ϵ-close distributions. Noise parameter γ represents depolarizing strength.

**Notes/Kinks:**
- Automatically inserts and removes measurement ops as needed.

---

## `xeb_truedist_empdist_noisy()`

**Purpose:**  
Computes XEB score between the ideal distribution and noisy empirical samples.

**Relation to Paper:**  
Used to measure Linear Cross Entropy Benchmark score between two distributions.

**Notes/Kinks:**
- Requires 10,000 to 1,000,000 shots for meaningful estimates (especially for n > 10).
- Negative values suggest sampling worse than uniform.

---

## `classical_fidelity(dist1, dist2)`

**Purpose:**  
Computes classical fidelity:  
![image](https://github.com/user-attachments/assets/09779258-ac1e-4407-b0d7-db4cb11e2e18)


**Relation to Paper:**  
Not directly used but related to ℓ₁-norm distance bounds.

**Notes/Kinks:**
- More sensitive to small differences than TVD.
- Requires positive values; zero-probability entries must be handled carefully.

---

## Special Considerations

### Bitstring Handling
- All functions use reversed (little-endian) bitstrings via `reverse_keys`.
- Qiskit’s native big-endian format would misalign with the paper’s conventions.

### Anti-Concentration Assumption
- Functions assume circuit depth ≥ Ω(log n) for valid metrics.
- Shallower circuits may not satisfy anti-concentration, leading to misleading results.

# Function Documentation  
*(Based on "A polynomial-time classical algorithm for noisy random circuit sampling")*

This documentation outlines utility functions used for analyzing noisy Random Circuit Sampling (RCS) simulations, aligned with statistical comparisons from the referenced paper. Metrics such as Total Variation Distance (TVD) and Cross Entropy Benchmarking (XEB) mirror tests used to compare classical and quantum output distributions.

---

### `compute_avg_xeb_varyingqubits()`

**Purpose:**  
Computes average XEB scores across various qubit counts using random circuits with depolarizing noise. Automates data generation for benchmarking claims from Theorem 1 (statistical indistinguishability).

**Relation to Paper:**  
Implements XEB tests from Section 1.2. Reflects the paper’s Appendix A which shows exponentially small XEB correlation in noisy RCS.

**Notes/Kinks:**  
- Requires `2^n` ideal probability calculations  
- Stores results as JSON to avoid recomputation  
- `γ` = per-gate depolarizing rate as defined in the paper

---

### `plot_avg_xeb_varyingqubits()`

**Purpose:**  
Visualizes how XEB decays with qubit count, contrasting ideal and noisy simulations.

**Relation to Paper:**  
Replicates Figure 1(e)-style results. Demonstrates `2^(-Θ(d))` XEB decay as predicted by Theorem 1.

**Notes/Kinks:**  
- Uses little-endian bitstrings via `reverse_keys`  
- Output filenames encode depth/noise settings for reproducibility

---

### `compute_avg_tvd()`

**Purpose:**  
Calculates average TVD between ideal and noisy output distributions across circuits.

**Relation to Paper:**  
Implements Eq. (26), demonstrating that noisy RCS distributions remain ε-close in TVD to classical outputs.

**Notes/Kinks:**  
- Supports both `log(n)` and linear depth  
- Parallelized over `gigaShots` for statistical reliability  
- Stores raw values for bootstrap-based error analysis

---

### `plot_avg_tvd()`

**Purpose:**  
Plots TVD scaling vs qubit count/depth, using bootstrapped error bars.

**Relation to Paper:**  
Reflects Section 3.1's `O(1)·e^(-2γℓ)` bound. Mimics visualizations like Figure 4.

**Notes/Kinks:**  
- Uses 3σ confidence from bootstrap samples  
- Compares real vs injected noise scenarios  
- Saves plots in publication-ready PDF/SVG

---

## Special Considerations

### **Verification Overhead**
- XEB needs `2^n` ideal calculations (classically exponential)  
- TVD becomes impractical for `n ≳ 20` due to exponential state space

### **Noise Simulation**
- Z-rotations emulate correlated/local phase noise  
- Uses γ = `1 - e^(-σ²/2)` as per the paper

### **Bitstring Handling**
- All functions use little-endian bitstrings via `reverse_keys`  
- Prevents mismatches with the paper’s theoretical expectations (unlike Qiskit’s big-endian format)

### References
Dorit Aharonov et al., "A polynomial-time classical algorithm for noisy random circuit sampling", arXiv:2211.03999
