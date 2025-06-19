This is a repo for implementing the Pauli Path Algorithm used for Noisy Random Circuit Sampling proposed in the [Aharonav et al.](https://arxiv.org/pdf/2211.03999) paper.

Our program uses three main components: Pauli Path generation, evaluation metrics, and Fourier coefficient calculation. It brings these pieces together in the class `NoisyProbDist`, which allows us to simulate a quantum experiment on our standard laptops.

### Table of Contents
- [Motivation](#motivation)
- [Overview of Research](#overview-of-research)
- [User Guides](#user-guides)
- [State of the Research](#state-of-the-research)

### Motivation
By leveraging the principles of quantum mechanics to process information, quantum computing has the potential to revolutionize fields such as cryptography, materials science, and drug discovery. Unfortunately, current quantum computers are not yet accurate enough to implement the groundbreaking algorithms developed for them. To unlock their potential, we need more research and investment in quantum computing.

To increase public support of quantum computing development, corporations often turn to a task known as random circuit sampling (RCS). Random circuit sampling can be efficiently carried out an a quantum computer, while it is generally believed that no classical computer (such as your laptop) could simulate it in polynomial time. As a result, companies such as Google tend to use it as a benchmark to establish quantum supremacy, which is the claim that quantum computers can outperform classical systems for specific tasks.

Recent work by Ahanorav et al. questions the supposed divide between classical and quantum computers for RCS. In their 2022 paper ["A polynomial-time classical algorithm for noisy random circuit sampling"](https://arxiv.org/pdf/2211.03999), the researchers present an efficient classical simulation for RCS. 

They share the theoretical algorithm and extend to future researchers the opportunity to determine how to practically implement the algorithm. If the algorithm can be efficiently implemented on a classical computer, then we must seek a task that is a more robust demonstration of quantum supremacy. We programmed the simultation outlined by [Ahanorav et al., 2022](https://arxiv.org/pdf/2211.03999) in python, implying the need of a new quantum supremacy benchmark.

## Overview of Research
This program classically simulates a quantum computing task known as noisy random circuit sampling. It accomplishes this by using an approach known as the Pauli path method. 

Before we outline the method, we define relevant words that will come up in the method description.
   >### Quantum Circuit Terminology
   >- **Qubit**: Quantum computers work with qubits, the counterpart of classical bits in quantum computing. Qubits exist in a superposition of classical states 0 and 1. This is mathematically represented as:
   ```math
   \big| \psi \rangle = \begin{bmatrix} \alpha \\ \beta \end{bmatrix}
   ```
   >where $\alpha,\beta \in \mathbb{C}$ are amplitudes determining the probabilities of measuring the qubit in the 0 state or 1 state, respectively. 
   >- **Born Rule**: The probability of measuring an outcome is the squared magnitude of its amplitude. In the above example, $|\alpha|^2$ is the liklihood of measuring 0 and $|\beta|^2$ is the probability of measuring 1.
   >- **Gate**: A quantum gate transforms qubits in such a way that they preserve the qubits' valid probability distribution (i.e. the probability of all outcomes summing to 1). Unlike a classical gate, a quantum gate must have the same number of outputs as there are inputs. Thus, the 2-qubit gates in our program take 2 qubits as input and output 2 qubits.
   >- **Unitary**: A matrix that mathematically specifies the operation of a quantum gate. To apply a gate to a qubit, we multiply its corresponding unitary with the qubit's vector. 
   >- **Measurement**: Measuring a qubit causes its superposition of 1 and 0 to collapse to exactly one of these outcomes. Measurement cannot be undone, so once you measure a qubit, there is no way to return it to its prior state.
   >- **Circuit**: Refers to a quantum circuit. Quantum circuits start with some fixed number of input qubits and then apply a sequence of gates on the input. The qubits may be measured at any point in the circuit, and upon measurement, collapse irreversibly to a particular outcome.
   >- **Random Circuit Sampling**: RCS is a benchmarking task designed to demonstrate quantum supremacy. The process involves repeatedly obtaining samples from the output distribution of randomly chosen quantum circuits. These circuits are characterized by an arbitrary set of gates and some fixed circuit architecture. In other words, the placement of the gates in the circuit is constant while the gates themselves are chosen arbitrarily.
   > - **Depth**: The number of gate layers in a quantum circuit. In the below example, it is 5.
   >- **Noise**: In the real world, quantum computers are susceptible to noise. Noise is errors caused by unexpected effects from neighboring qubits or external sources like radiation, magnetic field, electrical field, etc. As seen in the below figure, ideal RCS has no noise at all (a). Noise is indicated by the blue dots as seen in the noisy RCS diagram (b).
   > <p align="center">
   >   <img src="https://github.com/ilovelogic/Pauli-Path-Method/blob/main/images/RCS_circuit.png" width="800" />
   > </p>


   > ### Pauli Basis Terminology
   >- **Pauli**: The $2 \times 2$ Pauli matrices are denoted $I$, $X$, $Y$, and $Z$, and they comprise the Pauli basis.

   >- **Tensor Product**: A product defined in a way that, for unitaries $U \in \mathbb{C}^a$ and $V \in \mathbb{C}^b$, it preserves the property 
   ```math
   \left( U \otimes V\right)\left( v \otimes w\right)=\left( U v\right) \otimes \left( V w\right)
   ``` 
   >for all states $v \in \mathbb{C}^a$ and $w \in \mathbb{C}^b$.
   >- **Pauli operator**: For an $n$-qubit system, the corresponding Pauli operators are of the form
   ```math
   P = P_1 \otimes P_2 \otimes \cdots \otimes P_{n},
   ```
   >where $P_i \in \{ I,X,Y,Z\}$, for every $1 \leq i \leq n$. Put equivalently, Pauli operators are a tensor product of Pauli matrices. The prized property of Pauli operators is that every quantum state can be expressed as a linear combination of Pauli operators.
   >- **Hamming weight**: Counts the number of non-identity Paulis $\left( X,Y,Z\right)$ present in Pauli operators. For example, the Hamming weight of $I \otimes X \otimes Y$ is 2. When we discuss the Hamming weight of a Pauli path $s$, we denote it by $|s|$, and it is the sum of the number of non-identity Paulis in the Pauli operatrs of $s$.
   >- **Pauli path**: A Pauli path $s = (s_0, \ldots, s_d) \in \mathcal{P}_n^{d+1}$ is a sequence of Pauli operators representing the evolution of quantum states during a circuit's operation. Pauli paths are used to calculuate the Pauli path integral, which is a means of determining the output probability distribution of the circuit.
   

   With the basic terminology clarified, we define the Pauli path integral is defined as follows, according to the work of [Aharonav et al.](https://arxiv.org/pdf/2211.03999)

   >### Definition 1 (Pauli Path Integral)
   >Let $C = U_d U_{d-1} \cdots U_1$ be a quantum circuit acting on $n$ qubits, where $U_i$ is a layer of 2-qubit gates and $d$ is the circuit depth, and let $p(C, x)$ be the output probability distribution. The Pauli path integral is written as
   ```math
   p(C, x) = \sum_{s \in \mathbb{P}_n^{d+1}} f(C, s, x)
   ```
   >where $f(C,s,x)$ is the Fourier coefficient for output $x$ and Pauli path $s$, representing the contribution of that Pauli path in giving rise to outcome $x$.


Accordingly, to simulate a quantum circuit, we could construct all possible Pauli paths and use them to compute the above expression for all outcomes $x$ of our circuit. The result would be the theoretical probability distribution of our quantum circuit, as desired.

Determining all possible Pauli paths would be rather involved. Thankfully, there are some simplifications. To calculate the above expression, we only need to know all legal Pauli paths rather than all Pauli paths in general. Below, we give the definition of legal Paulis paths provided by [Aharonav et al.](https://arxiv.org/pdf/2211.03999).

   >### Definition 6 (Legal Pauli Path)
   >For a given circuit architecture, a Pauli path $s = (s_0, s_1, \ldots, s_d)$ is legal if the following two conditions are satisfied:
   >
   >1. For all 2-qubit gates in the circuit, its input and output Paulis are either both $II$, or both not $II$.
   >2. $s_0$ and $s_d$ only contain $I$ and $Z$.
   >
   > The reason for considering legal Pauli paths is that the illegal ones are irrelevant, as they contribute 0 to the Pauli path integral.

   Calculating all of the legal Pauli paths would still take an exponenetial amount of time, so it is impractical to take that approach. Instead we will capitalize on the noise present in current quantum computers.

   > ### Effect of noise
   > The amount that a Pauli path contributes to the output probability distribution of a noisy quantum circuit decreases exponentially as the Hamming weight of the Pauli path increases. Thus, we can estimate the output probability distribution of a noisy quantum circuit by only considering low weight paths.

Let $\gamma$ be the noise rate of our quantum circuit and $l$ be our upper bound on the Hamming weight of Pauli paths. Then we can approximate the probability distribution of outcome $x$ using

$$\tilde{p}(C, x) \approx \sum_{s \in \mathbb{P}_{n}^{d+1} : |s| \leq \ell} (1 - \gamma)^{|s|} f(C, s, x)$$

As a result, our code functions as follows:
- The `Path Generation` code generates all legal Pauli paths that fit the specified depth, number of qubits, gate positions, and upperbound on Hamming weight. It accomplishes this task using 5 classes, which are outlined in the `Path Generation` README. 
- Then the `Fourier Coefficient` code calculuates the Pauli path integral based off of the Pauli paths from the `Path Generation` code. 
- Lastly, `Brute Force RCS` verfies the resulting output distribution using statistical measurements XEB and TVD.

These componenets come together via the `NoisyProbDist` and `GetProbDist` classes, which are described in the first README file below.

## User Guides

### Legal Pauli Path Generation:
[Path Generation README](https://github.com/ilovelogic/Pauli-Path-Method/blob/main/Path_Generation/README.md)

### Fourier Coefficient Calculation:
[Fourier Coefficent README](https://github.com/ilovelogic/Pauli-Path-Method/blob/main/Pauli_Amplitude/README.md)

### Brute Force Simulation:
[Brute Force RCS README](https://github.com/ilovelogic/Pauli-Path-Method/blob/main/Brute_Force_RCS/README.md)

### Determining Probalities:
[Noisy Probability Distribution README](https://github.com/ilovelogic/Pauli-Path-Method/blob/main/Prob_Calc/README.md)

## State of the Research 

### Research Question:
Our original research question was "Is RCS a worthwhile approach to proving quantum supremacy? Should we focus our attempts elsewhere?" Our program significantly speeds up classical simulation runtime by the path truncation approach, so it appears that working with a different benchmark may be in order. If we want to demonstrate quantum supremacy, we need a task that is more resilient against classical simulation and enjoys a greater difference in runtime as compared to classcial computers.

Below, we outline what each one of us contributed to answering the research question.

### Pauli Path Generation (Anne Kelley):
I implemented legal Pauli path generation based on Lemma 8 of the [Aharonav et al.](https://arxiv.org/pdf/2211.03999) paper. Pauli path generation is essential to classically simulating RCS, given that you can't compute the probability distribution without tracing the Fourier coefficients of the paths. My code:
- Correctly generates all legal Pauli paths given the number of qubits, depth, circuit archicture, and upper bound on Hamming weight.
- Can encapsulate all possible paths in either tree format or list format. 
- Handles both 1D brickwork circuit architecture and 2D brickwork. 

Its main optimizations include:
- Grouping Pauli operators at each layer of the path generation so that you only needed to determine the next possible Pauli operators for one operator in each grouping, rather than for all of the Pauli operators at the layer.
- Building the tree version of all possible paths, which speeds up Fourier coefficient calculation.

A further line of research would be to generalize to the gate sets used by Google and USSTC, which are discussed in Section 4 of the [Aharonav et al.](https://arxiv.org/pdf/2211.03999) paper.

### Brute Force Simulation (Jesus Azpitarte):
Currently, random sampling is working for 1d brickwork circuits. Adopting 2D circuit generation for sampling shouldn't be difficult. It's merely making sure the labels work similarly as in the 1D case.

##### Explanation of 2D brickwork: 
The function builds a 2D brickwork quantum circuit on a rows × cols grid of qubits with a specified number of layers (depth). Each layer applies Haar-random 2-qubit gates in a staggered pattern: Types 1 and 2 apply gates horizontally on even and odd column pairs, while Types 3 and 4 apply them vertically on even and odd row pairs. These four types cycle every four layers, ensuring full connectivity over time. The resulting structure entangles qubits across both dimensions in a regular, non-overlapping pattern that supports universal quantum computation in a 2D architecture.

The evaluation metrics available are XEB, TVD, and Fideity (similar to TVD). Its recommended to use TVD as shallower circuits as the depth is not high enough to enable proper mixing/entangling of states. XEB is better for qubit ranges of aruond 10-20 qubits with log depth. XEB becomes more accurate with higher depth.

For the future, the brute force simulation should be adapted to 2D case.

### Fourier Coefficient and Marginal Sampling (Erika Lee):
I was responsible for implementing the Fourier coefficient calculation component of the project, which is the calculation computed for each Pauli path, which we use to produce the output distribution of noisy quantum circuits using the Pauli path framework. The Fourier coefficient module is now fully operational. It computes the quasi-probability amplitude f(C, s, x) for each legal Pauli path s = (s₀, ..., s_d), based on the expansion introduced in Aharonov et al. The implementation accounts for gate-local transitions by computing traces of the form Tr(sᵢ ⋅ Uᵢ ⋅ sᵢ₋₁ ⋅ Uᵢ†), and handles both full circuit output overlaps and partial marginal overlaps. Depolarizing noise is incorporated through a factor of (1 − γ)^|s|, where |s| denotes the Hamming weight of the path.

The code has been tested and validated on circuits of varying depth and architecture. It supports full compatibility with Pauli paths generated from the Lemma 8 module and works correctly on 1D architecture, and since it uses qubit indexing it should also work on 2D brickwork circuits though that has not been testsed yet. It also correctly aligns with Qiskit’s qubit ordering by reversing indices internally where needed. The Fourier coefficient routine is designed to be modular and composable, with functions for each subcomponent: input overlap, transition amplitudes per layer, and final output overlap.

As an extension of the fourier coefficient computation, there's also the task of marginal sampling. On the marginal sampling side, the functionality is partially implemented and needs to be completed in the future. It uses the Fourier framework to compute marginal probabilities over subsets of qubits by reordering summations, as described in Lemma 9 of the Aharonov paper. The sampling routine operates bit-by-bit, calculating q̄(C, x_i = 0 | fixed_bits) and q̄(C, x_i = 1 | fixed_bits) at each step and drawing the next bit from the resulting conditional distribution. This process allows us to achieve polynomial runtime sample full bitstrings in O(n) calls to the marginal routine instead of summing over all 2^n outcomes.