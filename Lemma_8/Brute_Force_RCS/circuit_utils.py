import math

import numpy as np
import matplotlib.pyplot as plt

from qiskit import QuantumCircuit, transpile, qpy
from qiskit.quantum_info import random_unitary, Statevector, Operator
from qiskit.visualization import plot_histogram

from qiskit_aer import AerSimulator
from qiskit_aer.primitives import SamplerV2 as Sampler
from qiskit_aer.noise import NoiseModel, depolarizing_error


def random_circuit(num_qubits:int, depth=None):
    """
    Produces a qiskit circuit with 1D brickwork architecture where each unitary
    is a Haar random unitary and with a measurement operator at the end.

    Args:
      num_qubits (int): The number of qubits to form the desired circuit. Must be > 3
      depth (int): The desired depth of the circuit. Defaults to log(n)

    Returns:
      qc: A Qiskit quantum circuit object with 1D brickwork architecture and a measurement.
    """
    qc = QuantumCircuit(num_qubits)
    
    if num_qubits < 3:
        print("Circuit must have 3 qubits or more")
        return qc

    # Specifying circuit depth is optional.
    if depth is None:
        depth = math.floor(math.log(num_qubits, 2))

    # Constructing interwieved layers of two qubit random unitaries
    for d in range(depth):
        if d % 2 == 0:  # Even layers 
            for i in range(0, num_qubits - 1, 2):
                random_gate = random_unitary(4).to_instruction()
                random_gate.label="random_unitary_layer_" + str(d)
                qc.append(random_gate, [i, i + 1])
        else:  # Odd layers
            for i in range(1, num_qubits - 1, 2):
                random_gate = random_unitary(4).to_instruction()
                random_gate.label="random_unitary_layer_" + str(d)
                qc.append(random_gate, [i, i + 1])

    return qc

def create_2d_brickwork_circuit(rows: int, cols: int, depth: int):
    """
    Create a 2D brickwork circuit with specified dimensions and depth.
    
    Args:
        rows (int): Number of rows in the 2D grid. Rows must be >= 2
        cols (int): Number of columns in the 2D grid
        depth (int): Circuit depth (number of layers)
    
    Returns:
        QuantumCircuit: A Qiskit quantum circuit with the 2d brickwork architecture where
        each unitary is a Harr random 2 qubit gate. 
        The circuit has a measurement operator at the end.
    """
    if (rows < 2):
        raise ValueError("Invalid input: rows can not be less than 2")
    
    num_qubits = rows * cols
    qc = QuantumCircuit(num_qubits)

    # Define layer types and their corresponding gate labels
    layer_types = [
        ("type1", "black"),
        ("type2", "red"),
        ("type3", "pink"),
        ("type4", "green")
    ]

    # Implement the brickwork pattern for the specified depth
    for layer in range(depth):
        current_layer_type, gate_label = layer_types[layer % 4]

        # Layer Type 1: 
        if current_layer_type == "type1":
            for r in range(0, rows, 1):  
                for c in range(0, cols - 1, 2):  # Staggered columns
                    q1 = r * cols + c
                    q2 = r * cols + c + 1
                    random_gate = random_unitary(4).to_instruction()
                    random_gate.label = gate_label
                    qc.append(random_gate, [q1, q2])


        # Layer Type 2:
        elif current_layer_type == "type2":
            for r in range(0, rows, 1):  
                for c in range(1, cols - 1, 2):  
                    q1 = r * cols + c
                    q2 = r * cols + c + 1
                    random_gate = random_unitary(4).to_instruction()
                    random_gate.label = gate_label
                    qc.append(random_gate, [q1, q2])

        # Layer Type 3:
        elif current_layer_type == "type3":
            for c in range(0, cols, 1):  
                for r in range(0, rows - 1, 2):  
                    q1 = r * cols + c
                    q2 = (r + 1) * cols + c
                    random_gate = random_unitary(4).to_instruction()
                    random_gate.label = gate_label
                    qc.append(random_gate, [q1, q2])

        # Layer Type 4:
        elif current_layer_type == "type4":
            for c in range(0, cols, 1):  
                for r in range(1, rows - 1, 2):  
                    q1 = r * cols + c
                    q2 = (r + 1) * cols + c
                    random_gate = random_unitary(4).to_instruction()
                    random_gate.label = gate_label
                    qc.append(random_gate, [q1, q2])

    return qc


def extract_gates_info(qc):
    """
    Extracts the position, matrix representation, and layer of each gate in a quantum circuit.
    
    Parameters:
    qc (QuantumCircuit): The quantum circuit to analyze.
    
    Returns:
    list of tuples: Each tuple contains (gate matrix, qubits applied to, layer number).
    """
    gate_info = []

    for instruction in qc.data:
        gate = instruction.operation
        qubits = instruction.qubits
        qubit_indices = [qc.find_bit(q).index for q in qubits]
        #gate_matrix = Operator(gate).data

        if gate.name != "measure":
            gate_matrix = Operator(gate).data


        # Parse layer number from label if possible
        label = gate.label
        if label and label.startswith("random_unitary_layer_"):
            layer = int(label.split("_")[-1])
        else:
            layer = None

        gate_info.append((gate_matrix, tuple(qubit_indices), layer))
    return gate_info

def run_ideal_simulation(qc: QuantumCircuit, shots=100):
    """
    Runs a quantum circuit using an ideal (noise-free) simulator and gathers
    outcome state 'shot' many times.

    Args:
        qc (QuantumCircuit) : Quantum circuit to simulate.
        shots (int): Number of times to sample the circuit execution.

    Returns:
        counts_ideal: A dictionary with measurement outcomes as keys and counts as values.
    """
    sampler = Sampler()
    job = sampler.run([qc], shots=shots)
    result_ideal = job.result()
    counts_ideal = result_ideal[0].data.meas.get_counts()

    return counts_ideal

def create_noise_model(depolarizing_param: float):
    """
    Creates a noise model with a depolarizing error applied to 2-qubit gates.
    Doesn't work for ideal depolarizing param of 0.

    Args:
        depolarizing_param (float): The depolarizing error rate. Must be between 0 and 1

    Returns:
        NoiseModel: A Qiskit noise model object with the specified noise.
    """

    if depolarizing_error > 1 or depolarizing_error <= 0:
        raise ValueError("depolarizing param should be between greater than 0 and less than or equal to 1")
    
    noise_model = NoiseModel()

    # Create a 2-qubit depolarizing error
    error = depolarizing_error(depolarizing_param, 2)

    # Attach the error to all gates labeled 'random_unitary'
    noise_model.add_all_qubit_quantum_error(error, 'random_unitary')

    # Add 'unitary' as a recognized basis gate in the simulation
    noise_model.add_basis_gates(['unitary'])

    return noise_model

# Gathers how many times a basis state is observed over 'shots' amount of times for a particular circuit
def run_noisy_simulation(circuit, noise_model, shots=100):
    sim_noise = AerSimulator(noise_model=noise_model)
    
    # Transpile the circuit for the noisy simulator
    tcirc_noise = transpile(circuit, sim_noise)
    
    # Run on the noisy simulator
    result = sim_noise.run(tcirc_noise, shots=shots).result()
    counts = result.get_counts(circuit)
    
    # Counts are helpful to return so we can visualize
    return counts

def run_noisy_simulation(circuit: QuantumCircuit, noise_model: NoiseModel, shots=100):
    """
    Executes a quantum circuit on a simulator with the specified noise model for
    'shots' amount of times.

    Args:
        circuit (QuantumCircuit): The quantum circuit to run.
        noise_model (NoiseModel): The noise model to apply.
        shots (int): Number of samples to collect from the circuit.

    Returns:
        counts: A dictionary with noisy measurement outcomes and their counts.
        measurement outcomes as keys and counts as values.
    """
    sim_noise = AerSimulator(noise_model=noise_model)

    # Transpile the circuit for execution on the noisy simulator
    tcirc_noise = transpile(circuit, sim_noise)

    # Run the simulation
    result = sim_noise.run(tcirc_noise, shots=shots).result()
    counts = result.get_counts(circuit)

    return counts

def complete_distribution(partial_dist: dict, num_qubits: int) -> dict:
    """
    Completes a partial probability distribution by filling in missing basis states with zero probability.

    Args:
        partial_dist (dict): A dictionary mapping binary basis states (as strings) to probabilities.
        num_qubits (int): The number of qubits in the system.

    Returns:
        dict: A complete distribution over all 2^n basis states, with missing states set to 0.0 probability.
    """
     
    num_states = 2 ** num_qubits
    bitstring_format = f'0{num_qubits}b'
    
    # Initialize full distribution with zero probability
    full_dist = {format(i, bitstring_format): np.float64(0.0) for i in range(num_states)}
    
    # Update with given probabilities
    for state, prob in partial_dist.items():
        full_dist[state] = np.float64(prob)
    
    return full_dist

def count_to_distribution(counts: dict[str, int], num_qubits, shots=100):
    """
    Converts raw counts into a probability distribution over all possible basis states.
    Basically for a outcome the probability is count/totalShots

    Args:
        counts (dict): Raw measurement counts from a circuit run.
        num_qubits (int): Number of qubits in the circuit.
        shots (int): Total number of circuit executions (samples).

    Returns:
        dict: A probability distribution over all 2^n basis states.
    """
    
    # Total number of basis states 
    num_states = 2 ** num_qubits

    # Format string to pad binary numbers (e.g., '03b' for 3 qubits → '001')
    bitstring_format = f'0{num_qubits}b'

    # Generate all basis states 
    basis_states = []
    for i in range(num_states):
        bitstring = format(i, bitstring_format)
        basis_states.append(bitstring)

    # Initialize distribution with zero probability for each state
    distribution = {state: np.float64(0.0) for state in basis_states}

    # Populate distribution with actual probabilities from counts
    for state, count in counts.items():
        distribution[state] = np.float64(count / shots)

    return distribution

def generate_emp_distribution(qc: QuantumCircuit, shots: int, noise=None, depth=None):
    """
    Generates the empirical probability distribution of a given quantum circuit 
    after simulating it under specified noise for shot amount of samples.

    Args:
        qc (QuantumCircuit): The circuit to simulate.
        shots (int): Number of measurement shots.
        noise (float, optional): Depolarizing noise parameter. Defaults to None (ideal).
        depth (int, optional): Not used in current function. Included for interface compatibility.

    Returns:
        dict: A probability distribution over bitstring outcomes.
    """
    
    # Run the simulation using either an ideal or noisy backend
    if (not noise or noise == 0.0):
        # Ideal simulation (no noise model applied)
        counts = run_ideal_simulation(qc, shots)
    else:
        # Noisy simulation using depolarizing noise model
        counts = run_noisy_simulation(qc, noise, shots)

    # Convert raw counts to empirical probability distribution
    distribution = count_to_distribution(counts, qc.num_qubits, shots)

    return distribution

