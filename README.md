This is a repo for implementing the Pauli Path Algorithm used for Noisy Random Circuit Sampling proposed in the Dorit Aharnov paper.

The project is separated into three main components, the Pauli Path generation, the regular brute force version + evaluation metrics, and the fourier coefficient calculation.

### Table of Contents
- [User Guides | READMES](#user-guides)
- [State of the Research](#state-of-the-research)

## User Guides
### Pauli Path Generation:
[Pauli Path Method - Lemma 8 README](https://github.com/ilovelogic/Pauli-Path-Method/tree/main/Lemma_8#readme)

### Brute Force Simulation:
[Brute Force RCS README](https://github.com/ilovelogic/Pauli-Path-Method/blob/main/Lemma_8/Brute_Force_RCS/README.md)

### Fourier Coefficient
[Fourier Coefficent README]()

## State of the Research 

#### Pauli Path:

### Brute Force Simulation State of the Research:
Currently, random sampling is working for 1d brickwork circuits. Adopting 2D circuit generation for sampling shouldn't be difficult. It's merely making sure the labels work similarly as in the 1D case.

##### Explanation of 2D brickwork: 
The function builds a 2D brickwork quantum circuit on a rows × cols grid of qubits with a specified number of layers (depth). Each layer applies Haar-random 2-qubit gates in a staggered pattern: Types 1 and 2 apply gates horizontally on even and odd column pairs, while Types 3 and 4 apply them vertically on even and odd row pairs. These four types cycle every four layers, ensuring full connectivity over time. The resulting structure entangles qubits across both dimensions in a regular, non-overlapping pattern that supports universal quantum computation in a 2D architecture.

The evaluation metrics available are XEB, TVD, and Fideity (similar to TVD). Its recommended to use TVD as shallower circuits as the depth is not high enough to enable proper mixing/entangling of states. XEB is better for qubit ranges of aruond 10-20 qubits with log depth. XEB becomes more accurate with higher depth.

For the future, the brute force simulation should be adapted to 2D case.

### Fourier Coefficient State of the Research


