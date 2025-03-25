# Lemma 8 Implementation: Toward a Polynomial-Time Classical Simulation of Noisy Random Circuit Sampling

A Python implementation of the algorithm described in Lemma 8 from the paper "[A polynomial-time classical algorithm for noisy random circuit sampling.](https://arxiv.org/pdf/2211.03999)" Constructs a list of all possible legal Pauli paths given the depth, number of qubits, gate positions, and upperbound on Hamming weight for a 2D architecture quantum circuit.

---

## Table of Contents

- [Introduction](#introduction)
- [Architecture Overview](#architecture-overview)
- [Class Documentation](#class-documentation)
  - [Circuit](#circuitpy)
  - [PauliPath](#pauli_pathpy)
  - [Layer](#layerpy)
  - [PauliOperator](#pauli_operatorpy)
  - [TestCircuit](#test_circuitpy)

---

## Introduction

The **Lemma 8 Implementation** plays a key role in classically simulating noisy random circuit sampling in polynomial time. 

---

## Architecture Overview

1. **Core Classes**:
   - `Circuit`: Constructs a list of all possible PauliPaths for a given circuit architecture and upperbound on Hamming weight.
   - `PauliPath`: Builds a list of Layer objects representing each layer of a circuit given a specific weight configuration.
   - `Layer`: Stores all of the PauliOperators at a specific layer in two hash maps, one with the lists sorted according to which PauliOperators propagate backward to the same list of PauliOperators and one sorted according to which PauliOperators propagate forward to the same PauliOperator list.
   - `PauliOperator`: Stores a Pauli operator and the list of all PauliOperators it can propagate forward to and a list of the ones to which it can propagate backward. Can generate the forward and/or backward list.
2. **Testing**:
   - `TestCircuit`: A suite of tests to validate functionality, both in general use cases and edge cases. Restricted to circuits with 0 to 25 qubits.

---

## Class Documentation

---

### circuit.py

**Overview**

**Initialization**
   '''Circuit(num_qubits:int, num_layers:int, l:int, gate_pos:List[List[tuple]])'''

**Attributes**
   - num_qubits
   - num_layers
   - gate_pos
   - max_weight
   - weight_combos
   - pauli_paths

**Methods**
   - init_pauli_paths()
   - enumerate_weights(weight_list:List[int], wiggle_room:int, num_layers_left:int)

---

### pauli_path.py

**Overview**

**Initialization**
   '''PauliPath(num_qubits:int, weight_combo:List[int],gate_pos:List[List[tuple]])'''

**Attributes**
   - num_qubits
   - depth
   - weight_combo
   - gate_pos
   - all_ops

**Methods**
   - build_min_configs()
   - unsorted_min_layer_ops(min_weight)
   - min_backward(min_layers,min_depth)
   - min_forward(min_layers,min_layer_ops,min_depth)
   - propagate_next(all_sibs:DefaultDict[tuple, List[PauliOperator]], pos_to_fill:DefaultDict[PauliOperator,List], backward:int, depth:int)

---

### layer.py

**Overview**

**Initialization**
   Layer(gate_pos:List[tuple]=None, backward:int=-1,pauli_ops:DefaultDict[tuple, List[PauliOperator]]=None)

**Attributes**
   - backward
   - gate_pos
   - forward_sibs
   - bacward_sibs
   - pos_to_fill
   - carry_over_qubits

**Methods**
   - check_qubits(unsorted_pauli_ops:List[PauliOperator])
   - find_sibs(unsorted_pauli_ops:List[PauliOperator])

---

### pauli_operator.py

**Overview**

**Initialization**
   PauliOperator(operator:List[str], backward_ops:List[PauliOperator] = None, forward_ops:List[PauliOperator] = None)

**Attributes**
   - operator
   - backward_ops
   - forward_ops
   - list_alloc

**Methods**
   - weight_to_operators(next_weight:int, pos_to_fill:List[tuple], backward:int)
   - list_allocs(num_p:int, num_w:int) = static
   - append_to_layers(sibs:List[PauliOperator],indices:tuple, strs:tuple, r_start:int, r_end:int) = static
   - add_gate_input(sibs:List[PauliOperator], num_RRs:int, pos_to_fill:List[tuple], r_start:int)

---

## test_circuit.py
