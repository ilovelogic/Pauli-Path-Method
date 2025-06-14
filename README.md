This is a repo for implementing the Pauli Path Algorithm used for Noisy Random Circuit Sampling proposed in the Dorit Aharnov paper.

The project is separated into three main components, the Pauli Path generation, the regular brute force version + evaluation metrics, and the fourier coefficient calculation.

### Table of Contents
- [Overview of Research](#overview-of-research)
- [User Guides](#user-guides)
- [State of the Research](#state-of-the-research)

## Overview of Research
This program plays a key role in classically simulating noisy random circuit sampling in polynomial time. It accomplishes this by using an approach known as the Pauli path method. 

Before we outline the method, we define relevant words that will come up in the method description.
   >### Quantum Circuit Terminology
   >- **Qubit**: Quantum computers work with qubits, the counterpart of classical bits in quantum computing. Qubits exist in a superposition of classical states 0 and 1. This is mathematically represented as: $$\ket{\psi} = \alpha\ket{0} + \beta\ket{1} = \begin{bmatrix} \alpha \\ \beta \end{bmatrix}$$, where $\alpha$ and $\beta$ are complex amplitudes determining the probabilities of measuring the qubit in $\ket{0}$ or $\ket{1}$.
   >- **Gate**: In the context of our program, gate refers to a quantum gate. A quantum gate transforms qubits in such a way that they preserve the qubits' valid probability distribution (i.e. the probability of all outcomes summing to 1). Unlike a classical gate, a quantum gate must have the same number of outputs as there are inputs. Thus, the 2-qubit gates in our program take 2 qubits as input and output 2 qubits.
   >- **Unitary**: A matrix that mathematically specifies the operation of a quantum gate.
   >- **Measurement**: Measuring a qubit causes its superposition of 1 and 0 to collapse to exactly one of these outcomes. Measurement cannot be undone, so once you measure a qubit, there is no way to return it to its prior state.
   >- **Circuit**: Refers to a quantum circuit. Quantum circuits start with some fixed number of input qubits and then apply a sequence of gates on the input. The qubits may be measured at any point in the circuit, and upon measurement, collapse irreversibly to a particular outcome.
   >- **Random Circuit Sampling**: RCS is a benchmarking task designed to demonstrate quantum supremacy. The process involves repeatedly obtaining samples from the output distribution of randomly chosen quantum circuits. These circuits are characterized by an arbitrary set of gates and some fixed circuit architecture. In other words, the placement of the gates in the circuit is constant while the gates themselves are chosen arbitrarily.
   >- **Noise**: In the real world, quantum computers are susceptible to noise. Noise is errors caused by unexpected effects from neighboring qubits or external sources like radiation, magnetic field, electrical field, etc. As seen in the below figure, ideal RCS has no noise at all (a). Noise is indicated by the blue dots as seen in the noisy RCS diagram (b).
   >- **Depth**: The number of gate layers in a quantum circuit. In the above circuit example, the depth is 5.


   > ### Pauli Basis Terminology
   >- **Pauli**: The four Pauli matrices comprise the Pauli basis. They are as follows:
   >$$I = 
   >\begin{pmatrix}
   >1 & 0 \\
   >0 & 1
   >\end{pmatrix}, \quad
   >X = 
   >\begin{pmatrix}
   >0 & 1 \\
   >1 & 0
   >\end{pmatrix}, \quad
   Y = 
   \begin{pmatrix}
   0 & -i \\
   i & 0
   \end{pmatrix}, \quad
   Z = 
   \begin{pmatrix}
   1 & 0 \\
   0 & -1
   \end{pmatrix}.$$
   >- **Tensor Product**: A product defined in such a way that, for matrices $U \in \mathbb{C}^a$ and >$V \in \mathbb{C}^b$, it preserves the property 
   >$$\left( U \otimes V\right)\left( v \otimes w\right)=\left( U v\right) \otimes \left( V w\right)$$ 
   >for all states $v \in \mathbb{C}^a$ and $w \in \mathbb{C}^b$.
   >- **Pauli operator**: For an $n$-qubit system, the corresponding Pauli operators are of the form
   >$$P=P_1 \otimes P_2 \otimes \dots \otimes P_n$$
   >where $P_i \in \{ I,X,Y,Z\}$, for every $1 \leq i \leq n$. Put equivalently, Pauli operators are a tensor product of Pauli matrices. The prized property of Pauli operators is that every quantum state can be expressed as a linear combination of Pauli operators.
   >- **Hamming weight**: Counts the number of non-identity Paulis ($X,Y,Z$) present in Pauli operators. For example, the Hamming weight of $I \otimes X \otimes Y$ is 2. When we discuss the Hamming weight of a Pauli path $s$, we denote it by $|s|$, and it is the sum of the number of non-identity Paulis in the Pauli operatrs of $s$.
   >- **Pauli path**: A Pauli path $s = (s_0, \ldots, s_d) \in \mathcal{P}_n^{d+1}$ is a sequence of Pauli operators representing the evolution of quantum states during a circuit's operation. Pauli paths are used to calculuate the Pauli path integral, which is a means of determining the output probability distribution of the circuit.
   

   With the basic terminology clarified, we define the Pauli path integral is defined as follows, according to the work of [Aharonav et al.](https://arxiv.org/pdf/2211.03999)

   >### Definition 1 (Pauli Path Integral)
   > Let $C = U_d U_{d-1} \cdots U_1$ be a quantum circuit acting on $n$ qubits, where $U_i$ is a layer of 2-qubit gates and $d$ is the circuit depth. The Pauli path integral is written as: 
   >
   > $$p(C, x) = \sum_{s_0, \ldots, s_d \in \mathcal{P}_n} \text{Tr}(|x\rangle \langle x| s_d) \, \text{Tr}(s_d U_d s_{d-1} U_d^\dagger) \cdots \text{Tr}(s_1 U_1 s_0 U_1^\dagger) \, \text{Tr}(s_0 | 0^n \rangle \langle 0^n |).$$

Note that $p(C, x)$ is the output probability distribution for outcome $x$. Accordingly, to simualte a quantum circuit, we could construct all possible Pauli paths and use them to compute the above expression for all outcomes $x$ of our circuit.

Determining all possible Pauli paths seems rather involved. Thankfully, there are some simplications. To calculate the above expression, we only need to know all legal Pauli paths rather than all Pauli paths in general. Below, we define a legal Pauli path.

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
- The `Lemma 8` code generates all legal Pauli paths that fit the specified depth, number of qubits, gate positions, and upperbound on Hamming weight. It accomplishes this task using 5 classes, which are outlined in the `Lemma 8` README. 
- Then the `Fourier Coefficient` code calculuates the Pauli path integral based off of the Pauli paths from the `Lemma 8` code. 
- Lastly, `Brute Force RCS` verfies the resulting output distribution using statistical measurements XEB and TVD.

These componenets come together via the `NoisyProbDist` and `GetProbDist` classes, which are described in the first README file below.

## User Guides

### Determining Probalities:
[Noisy Probability Distribution README](https://github.com/ilovelogic/Pauli-Path-Method/tree/main/Lemma_8/Pauli_Amplitude/ProbCalcREADME/)

### Pauli Path Generation:
[Lemma 8 README](https://github.com/ilovelogic/Pauli-Path-Method/tree/main/Lemma_8#readme)

### Fourier Coefficient Calculation:
[Fourier Coefficent README]()

### Brute Force Simulation:
[Brute Force RCS README](https://github.com/ilovelogic/Pauli-Path-Method/blob/main/Lemma_8/Brute_Force_RCS/README.md)
\
\
## State of the Research 

### Research Question:
Our original research question was "Is RCS a worthwhile approach to proving quantum supremacy? Should we focus our attempts elsewhere?" Our program significantly speeds up classical simulation runtime by the path truncation approach, so it appears that working with a different benchmark may be in order. If we want to demonstrate quantum supremacy, we need a task that is more resilient against classical simulation and enjoys a greater difference in runtime as compared to classcial computers.

Below, we outline what each one of us contributed to answering the research question.

### Anne Kelley:
I handled generating the legal Pauli paths restricted by an upper bound on Hamming weight. The Pauli path generation is essential to classically simulating RCS, given that you can't compute the probability distribution without tracing the Fourier coefficients of the paths. My code:
- Correctly generates all legal Pauli paths given the number of qubits, depth, circuit archicture, and upper bound on Hamming weight.
- Can encapsulate all possible paths in either tree format, which speeds up Fourier coefficient calculations, or list format. 
- Handles both 1D brickwork circuit architecture and 2D brickwork. 

A further optimization would be to generalize to the gate sets used by Google and USTC, which are discussed in Section 4 of the [Aharonav et al.](https://arxiv.org/pdf/2211.03999) paper.

### Brute Force Simulation State of the Research (Jesus Azpitarte):
Currently, random sampling is working for 1d brickwork circuits. Adopting 2D circuit generation for sampling shouldn't be difficult. It's merely making sure the labels work similarly as in the 1D case.

##### Explanation of 2D brickwork: 
The function builds a 2D brickwork quantum circuit on a rows × cols grid of qubits with a specified number of layers (depth). Each layer applies Haar-random 2-qubit gates in a staggered pattern: Types 1 and 2 apply gates horizontally on even and odd column pairs, while Types 3 and 4 apply them vertically on even and odd row pairs. These four types cycle every four layers, ensuring full connectivity over time. The resulting structure entangles qubits across both dimensions in a regular, non-overlapping pattern that supports universal quantum computation in a 2D architecture.

The evaluation metrics available are XEB, TVD, and Fideity (similar to TVD). Its recommended to use TVD as shallower circuits as the depth is not high enough to enable proper mixing/entangling of states. XEB is better for qubit ranges of aruond 10-20 qubits with log depth. XEB becomes more accurate with higher depth.

For the future, the brute force simulation should be adapted to 2D case.

### Fourier Coefficient State of the Research


