# Pauli Path Generation

A Python implementation of the algorithm described in Lemma 8 from the work of [Aharonav et al.](https://arxiv.org/pdf/2211.03999). Constructs a list of all possible legal Pauli paths given the depth, number of qubits, gate positions, and upperbound on Hamming weight for a 2D architecture quantum circuit.

---

## Table of Contents

- [Introduction](#introduction)
- [Architecture Overview](#architecture-overview)
- [Class Documentation](#class-documentation)
  - [PauliOperator](#paulioperator)
  - [PauliOpLayer](#paulioplayer)
  - [PauliPathTrav](#paulipathtrav)
  - [XYZGeneration](#xyzgeneration)
  - [CircuiSim](#circuitsim)

---

## Introduction

The `Lemma_8` program generates all legal Pauli paths that fit the specified depth, number of qubits, gate positions, and upperbound on Hamming weight. It accomplishes this task using 5 classes, which are outlined below.

---

## Architecture Overview

1. **Core Classes**:
   - [PauliOperator](#paulioperator): Stores a Pauli operator, the list of all `PauliOperator` objects that can precede this operator in a legal Pauli path, and a list of those that can come after it. Can generate both of these lists. Each Pauli operator is represented as a list of strs, where each str represents which Pauli matrices can be at that particular index in the tensor product comprising the Pauli operator.
   - [PauliOpLayer](#paulioplayer): Keeps track of all the `PauliOperator` objects that can be the ith Pauli operator in a legal Pauli path, restricted by the circuit architecture and weight configuration. Uses two hash maps, one containing lists sorted according to which `PauliOperator` objects propagate backward to the same list of `PauliOperator` objects and one sorted according to which `PauliOperator` objects propagate forward to the same `PauliOperator` list.
   - [PauliPathTrav](#paulipathtrav): For traversing different possibile branches of our Pauli path. Builds a list of `PauliOpLayer` objects, where the ith `PauliOpLayer` in the list contains all the possibilities for the ith Pauli operator of the Pauli path.
   - [CircuitSim](#circuitsim): Constructs a list of all possible `PauliPathTrav` objects for a given circuit architecture and upperbound on Hamming weight.
   - [XYZGeneration](#xyzgeneration): Builds a list (`next_gen`) of all `XYZGeneration` objects that can come after this `XYZGeneration` to form valid Pauli path traversals.

---

## Class Documentation

---

### PauliOperator

**List of Strs Representation**\
Recall that each Pauli operator is a tensor product of matrices drawn from the $2 \times 2$ Paulis $X$, $Y$, $Z$, and $I$. Accordingly, the `PauliOperator` class represents a Pauli operator by a list of strs (the `operator` attribute), where the ith str in the list characterizes the ith Pauli in the tensor product. 

When we first propogate from a `PauliOperator` object, we use the following strs in `operator` attributes. 
   - "I": The associated Pauli matrix is the $2 \times 2$ identity matrix.
   - "R": The Pauli matrix can be either $X$, $Y$, or $Z$ (any of the non-identity Paulis).
   - "N": The Pauli matrix must be the same as the Pauli at the same index in the next Pauli operator in the path. We use "N" in the case where there is a non-identity Pauli whose associated qubit is not sent into any gates between this operator and the next, so it must not change between layers.
   - "P": The Pauli matrix is specified by whatever the Pauli is at the same index in the prior Pauli operator in the path. This is for the case of a non-identity Pauli whose related qubit does not enter any gates between this Pauli operator and the prior one.

We use "R", "N", and "P" rather than directly including "X", "Y", and "Z" in order to avoid blowing up in memory at this stage. Since each "R" encapsulates the possibility of "X", "Y", or "Z", our number of Pauli paths to store is exponentially less than it would be if we opted 
to store each individual combination of choosing either "X", "Y", or "Z" for all of "R"s in our list.

Later, using the class `XYZGeneration`, we will replace each `PauliOperator`'s "R", "N", and "P" with "X", "Y", and "Z", yet we will do so in a tree-like structure that handles the memory problem. \
\
<img src="/workspaces/Lemma-8/images/PauliOperator_ops.png" width="450" />\
**Immediate Neighbors**\
It is also essential that we situate each Pauli operator in terms of its neighbors in a legal Pauli path. Thus, we store a reference (`prior_ops`) to the `PauliOperator` objects that can directly precede this `PauliOperator` in a legal Pauli path. Likewise, we maintain an attribute (`next_ops`) for the `PauliOperator` objects that can come directly after this `PauliOperator`.

**Initialization**\
   `PauliOperator(operator:List[str], prior_ops:List[PauliOperator] = None, next_ops:List[PauliOperator] = None):`
   >Instantiates a PauliOperator object with its `operator` attribute initialized, and its `prior_ops` and `next_ops` attributes initialized if those parameters were included in the call.

**Attributes**
   - `operator`: A list of strs representing this `PauliOperator`, eiither in terms of "I"s, "R"s, "N"s, and "P"s or "I"s, "X"s, "Y"s, and "Z"s.
   - `prior_ops`: A list containing all `PauliOperator` objects that can directly precede this `PauliOperator` in a legal Pauli path.
   - `next_ops`: A list containing each of the `PauliOperator` objects that could come directly after this `PauliOperator` in a legal Pauli path.
   - `list_alloc`: A 2D array, where `list_alloc[i,j]` is the number of ways we
   can fill i gates with non-identity I/O using an overall Hamming weight of j

**Methods:**

- **`weight_to_operators(sib_ops: List[PauliOperator], next_weight: int, pos_to_fill: List[tuple], backward: int)`**  
  If `backward` is 1, this method determines all possible `PauliOperator` objects that can directly precede the given `PauliOperator` in a legal Pauli path. These objects are appended to the `prior_ops` attribute of the class, with their `operator`s being in terms of "I", "R", "N", and "P". Otherwise, it determines all possible `PauliOperator` objects that can directly follow the given one in a legal Pauli path and appends them to the `next_ops` attribute.

- **`list_allocs(num_p: int, num_w: int)`**  
  A static method that calculates the number of ways to distribute `num_w` (the total Hamming weight) across `num_p` non-identity gate positions. It returns a 2D integer array where each entry `(i, j)` represents the number of ways to fill `i` gates with a total Hamming weight of `j`.

- **`find_next_operators(sibs: List[PauliOperator], num_RRs: int, pos_to_fill: List[tuple], r_start: int)`**  
  A helper method of `weight_to_operators` which recursively fills the next non-identity I/O gate positions with either "I" and "R", "R" and "I", or "R" and "R", until we reach the bases case where `num_RRs` is 0 or the number of positions left to fill is equal to `num_RRs`.

- **`edit_ops(sibs: List[PauliOperator], indices: tuple, strs: tuple, r_start: int, r_end: int)`**  
  A static method that fills the first index of `indices` with the first str of `strs` and the second index with the second str for the `operator` attribute of all `PauliOperator` objects of `sibs` in the range [`r_start`,`r_end`).

---

### PauliOpLayer

<img src="/workspaces/Lemma-8/images/backward_forward_sibs.png" width="400" />

**Sorting for Propagation**\
We determine all possibile Pauli paths by working terms of layers, which are represented by the object `PauliOpLayer`. Each `PauliOpLayer` stores all the `PauliOperator` objects which could be the Pauli operator for some particular position in a legal Pauli path. 

We group `PauliOperator` objects at this `PauliOperator` layer in lists according to which propagate backward to the same list of `PauliOperator` objects and which propagate forward together.

We store each backward-twinning list in a DefaultDict (`backward_rnp_sibs`), where the list is the value and the key is a tuple of the list of non-identity gate positions between any `PauliOperator` in the list and all of the `PauliOperator` objects in the list it propagates backward to, along with a List of strs representing the non-gate qubits of the `PauliOperator` in this list in the context of looking backward. `forward_rnp_sibs` is a DefaultDict with the same set up except that it stores lists containing `PauliOperator` objects that propagate forward to the same list of possibilities at the *next* layer.

**Initialization**\
   `PauliOpLayer(gate_pos:List[tuple]=None, backward:int=-1,pauli_ops:DefaultDict[tuple, List[PauliOperator]]=None)`
   >Takes as argument `pauli_ops`, a `DefaultDict` that is either sorted according to backward or forward propagation at this layer of the Pauli path. Calls methods `check_qubits` and `group_sibs` to group the `PauliOperator` objects in `pauli_ops` based off of the opposite propagation direction. Assigns attribute `forward_rnp_sibs` to be the `DefaultDict` sorted by which propagate forward together and `backward_rnp_sibs` to be the `DefaultDict` that groups based on backward propagation.

**Attributes**
   - `backward`: An int that is 1 if we need to propagate backward from this `PauliOpLayer` and 0 otherwise.
   - `gate_pos`: A list of int tuples, where each int tuple stores the two I/O indices of a gate between this PauliOpLayer and either the prior PauliOpLayer (if backward) or the next `PauliOpLayer` (otherwise).
   - `pos_to_fill`: A DefaultDict whose keys are the `PauliOperator` objects at this `PauliOpLayer` and whose values are lists of int tuples containing the non-identity I/O gate positions between the `PauliOperator` key and the `PauliOpLayer` to which we're propagating.
   - `forward_rnp_sibs`: A DefaultDict that sorts all the `PauliOperator` objects at this `PauliOpLayer` according to their having matching gate positions with non-identity I/O and the same non-gate qubits when propagating forward. The non-identity I/O gate positions list and the list of non-gate qubit strs are both converted into tuples and used as the key for the DefaultDict. The list of "family" PauliOperator objects is the value for the DefaultDict.
   - `bacward_sibs`: A DefaultDict with the same setup as `forward_rnp_sibs` except that it is sorted according to matching the attributes when propagating backward. Both `forward_rnp_sibs` and `backward_rnp_sibs` enjoy the property that each of their "family" lists contain a grouping of `PauliOperator` objects that all propagate to the same list of `PauliOperator` objects at a neighboring circuit Layer.
   - `carry_over_qubits`: A list of lists of strs, where each list of strs is a copy of one the `PauliOperator` object's operator at this `PauliOpLayer` with the edit that its qubits at gate positions are all set to "I". This allows us to check equality of non-gate qubits for `PauliOperator` objects by simply comparing their lists in `carry_over_qubits`. 

**Methods**
   - **`check_qubits(unsorted_pauli_ops:List[PauliOperator])`**\
   Uses the attribute `gate_pos` to determine which gate positions have non-identity I/O, for each `PauliOperator` in `unsorted_pauli_ops`, in order to fill out the attribute `pos_to_fill`. Also uses `gate_pos` to correctly edit `carry_over_qubits` so that each of its entries only characterizes non-gate qubits, for each `PauliOperator` in `unsorted_pauli_ops`. 

   - **`group_sibs(unsorted_pauli_ops:List[PauliOperator])`**\
   Relies on the information obtained from calling `check_qubits` to sort all the `PauliOperator` objects of this `PauliOpLayer` into lists according to which have the same set of non-gate qubits and non-identity I/O gate positions. Accomplishes that using a DefaultDict where the keys are tuples comprised of a `PauliOperator` object's associated entry of `pos_to_fill` and `carry_over_qubits`. 

---

### PauliPathTrav

**Overview**\
<img src="/workspaces/Lemma-8/images/PauliPath_layers.png" width="800" />
Given a circuit architecture and weight configuration, we generate a list of `PauliOpLayer` objects (`layers`), which encapsulates all the possibilities of what `PauliOperator` objects we are free to use at each index of a legal Pauli path. 

**Initialization**\
   `PauliPathTrav(num_qubits:int, weight_combo:List[int],gate_pos:List[List[tuple]])`
   >Begins by creating a `PauliOpLayer` that stores all possible `PauliOperator` objects meeting the weight requirement at the minimum weight layer. Then propagates backward from the min layer until reaching the beginning of the path and forward from the min layer until getting to the end of the path. As it propagates, it builds the `PauliOpLayer` for each index of the Pauli path and stores it in the corresponding index of the attribute `layers`. 

**Attributes**
   - `num_qubits`: An int that is the number of qubits in the circuit.
   - `num_op_layers`: An int representing the number of Pauli operators in any valid Pauli path represented by this `PauliPathTrav`
   - `weight_combo`: A list of ints, where the ith list entry is the required Hamming weight for any ith Pauli operator in a Pauli path of this `PauliPathTrav`.
   - `gate_pos`: A list of lists of int tuples, where the ith list of int tuples contains the gate positions between the ith and i+1st Pauli operators.
   - `layers`: A list of `PauliOpLayer` objects, where the ith `PauliOpLayer` object keeps track of all possible ith Pauli operators in our Pauli path.

**Methods**
   - **`build_min_configs()`**\
   Determines and returns the index of the lowest Hamming weight layer and all possible `PauliOperator` objects at that layer. Also returns the non-identity gate positions between the min weight layer and its prior layer and the non-identity gate positions between the min layer and its next layer, both as `List[tuple]`. Uses helper methods `unsorted_min_layer_ops`, `min_backward`, and `min_forward`.

   - **`unsorted_min_layer_ops(min_weight)`**\
    Returns an unsorted list of all PauliOperators whose number of qubits match the attribute `num_qubits` and whose Hamming weight equals `min_weight`.

   - **`min_backward(min_layers:List[PauliOperator],min_depth:int)`**\
   Instantiates the `PauliOpLayer` object `min_layer_ops` and sets its `backward_rnp_sibs` attribute, which is a `DefaultDict` whose values are `PauliOperator` objects (in "R", "N", "P, and "I") that are grouped by the key of the tuple of their gate positions and non-gate qubits in the backward direction. As a result, `backward_rnp_sibs` stores the `PauliOperator` objects at index `min_index` in the Pauli path sorted according to which propagate *backward* to the same list of `PauliOperator` objects. Returns `backward_rnp_sibs` and a `DefaultDict` that, for every `PauliOperator` possible at this layer, keeps track of all positions in the `PauliOperator` immediately prior to this `PauliOperator` that will require non-identity Paulis.

   - **`min_forward(min_layers:List[PauliOperator],min_layer_ops:PauliOplayer,min_depth:int)`**\
   Given the `PauliOpLayer` object `min_layer_ops`, instantiates its `forward_rnp_sibs` attribute. Its `forward_rnp_sibs` is a `DefaultDict` with values of `PauliOperator` objects categorized by their gate positions and non-gate qubits in the context of moving forward in the circuit. As a result, `forward_rnp_sibs` stores the `PauliOperator` objects at index `min_index` in the Pauli path organized by which propagate *forward* to the same list of `PauliOperator` objects. Returns `forward_rnp_sibs` and the `DefaultDict` that stores, for every `PauliOperator` at this layer, all qubit indices in the `PauliOperator` directly preceding this `PauliOperator` that need non-identity Paulis.

   - **`propagate_next(all_sibs:DefaultDict[tuple, List[PauliOperator]], pos_to_fill:DefaultDict[PauliOperator,List], backward:int, depth:int)`**\
   Takes in a list of forward or backward sibling operators at a layer (`all_sibs`), and determines the new sibling operators that each sibling operators list in the input list propagate to. Uses helper function weight_to_operaters from the PauliOperator class to get the sibling operators that all the PauliOperators in any given sibling operators of the input propagate to.
   

---

### XYZGeneration

<img src="/workspaces/Lemma-8/images/xyz_gen.png" width="700" />\

**Nesting Structure**\
The `XYZGeneration` class uses a recursive structure to generate Pauli paths. The constructor takes as input a list of `PauliOperator` objects (`parent_ops`), which all share the same list of `PauliOperator` objects that could come after them in a legal Pauli path. It also takes a `List[PauliOperator]` (`pauli_path`), which stores the current Pauli path, and the parameter `next_index` lets us know which index will hold the `PauliOperator` object that directly comes after one of the Pauli operators of the `parent_ops` list. 

**Initialization**\
   `XYZGeneration(pauli_ops:List[PauliOperator],next_index:int,pauli_path:List[PauliOperator])`

**Attributes**
   - `parent_ops`: The list of `PauliOperator` objects for a particular index in the Pauli path that have the same selection from "X", "Y", and "Z" for their non-gate non-identity qubits.
   - `next_gen`: A list of `XYZGeneration` objects, where the `parent_ops` attribute of each of these `XYZGeneration` contains all the `PauliOperator` objects that could come directly after any of the `PauliOperator` objects in a valid Pauli path.

**Methods**
   - **`rnp_to_xyz(next_index:int, pauli_path:List[PauliOperator])`**\


   - **`fill_pos_lists(next_op:PauliOperator, r_pos_list: List[int], n_pos_list: List[int])`**\
   Fills `r_pos_list` with all "R" qubit positions and `n_pos_list` with all non-carry "N" qubit positions in `next_op`. Note that an "N" that carries is a qubit position that remains a non-gate position in every following `PauliOperator` in the Pauli path. Replaces any "P"s it encounters with the non-identity Pauli at the same index in the immediately preceding `PauliOperator`.

   - **`fill_in_pos(filled_pos:List[PauliOperator],pos_list:List[int], pauli:str, index:int, start:int)`**\
   Updates `filled_pos` to contain all `PauliOperator` objects whose qubits at the positions in specified by `pos_list` hold all possible combinations of "X", "Y", and "Z".
   
   - **`carries_to_the_end(pauli_path_index:int, i:int)`**\
   Checks if the non-gate qubit at index `i` in the Pauli operator at index `pauli_path_index`
   in the Pauli path carries to the end of the Pauli path. In other words, investigates if it remains a non-gate qubit until the last layer of the Pauli path. In this is the case, that position in that Pauli operator and onward in the path would be forced to be a "Z" to meet the all "I"s and "Z" at last layer requirement.

   - **`rp_to_z(next_op:PauliOperator, pauli_path:List[PauliOperator])`**\
   Replaces all "R"s and "P"s in the last layer with "Z"s to meet the restriction of only "I"s and "Z" being present in the last `PauliOperator` of the Pauli path.
  

---

### CircuitSim

**Overview**


The `CircuitSim` class represents a classical simulation of a noisy random circuit. Its constructor initiailizes a list of all possible `PauliPathTrav` objects, given the circuit architecture and an upperbound on Hamming weight.

**Initialization**\
   `CircuitSim(num_qubits:int, l:int, gate_pos:List[List[tuple]])`

**Attributes**
   - `num_qubits`: An int that is the number of qubits in the circuit.
   - `num_op_layers`: An int that is the number of Pauli operators in any of the circuit's valid Pauli paths, which equals the number of gates in the circuit plus 1.
   - `gate_pos`: A list of list of int tuples, where the ith list of int tuples represents all the gate positions in the ith gate layer of the circuit.
   - `max_weight`: An int that is the upper bound on the total Hamming weight of any Pauli path used for the simulation.
   - `weight_combos`: A list of lists of ints, where each list of ints represents an indexed assignment of weights to Pauli operators in a legal Pauli path.
   - `pauli_paths`: A list of all possible PauliPathTrav objects, given the upper bound on Hamming weight and circuit architecture.

**Methods**
   - **`valid_gate_pos(num_qubits:int, gate_pos:List[List[tuple]])`**\
    A static method that checks whether the specified number of qubits (`num_qubits`) and gate position array (`gate_pos`) comprise a valid circuit architecture.

   - **`enumerate_weights(weight_list:List[int], wiggle_room:int, num_layers_left:int)`**\
    Recursively fills the `weight_combo` attribute of the `CircuitSim` with each distinct list that specifies the Hamming weights of a legal Pauli path, taking into account the restrictions of the circuit architecture. For each list, the ith int in the list assigns the Hamming weight of the ith Pauli operator of the Pauli path.

   - **`init_pauli_paths()`**\
    Initiates the list of all `PauliPathTrav` objects that fit the circuit architecture and upper bound on Hamming weight.

   - **`travs_to_list()`**\
   Translates each `PauliPathTrav` in `pauli_path_travs` into a list of Pauli paths, still in "R", "N", "P", and "I", and stores these lists in attribute `rnp_pauli_paths`. Calls helper function `trav_to_list(trav:PauliPathTrav)` on each `PauliPathTrav` object to obtain its corresponding list of paths.

   - **`trav_to_list(trav:PauliPathTrav)`**\
   For every `PauliOperator` in the `forward_rnp_sibs` attribute of `trav`, calls on helper function `pauli_op_hopping`. Uses `pauli_op_hopping` to append all Pauli paths starting with that `PauliOperator` and taking some traversal through `trav` to the `List[List[PauliOperator]]` `trav_list`. Returns `trav_list` once it is filled.

   - **`trav_branching(trav_list:List[List[PauliOperator]], partial_pauli_path:List[PauliOperator], pauli_op:PauliOperator)`**\
   Recursively traverses every Pauli path through a `PauliPathTrav` and appends them to `trav_list`. The initial call must be on an empty `partial_pauli_path` and one of the `PauliOperator` objects in the first layer of a `PauliPathTrav`.

   - **`build_xyz_trees()`**\
   For each Pauli path list in `rnp_pauli_paths`, constructs its corresponding `XYZGeneration` tree, with the tree's branching representing different valid selections of"X", "Y", and "Z". Stores the root of each tree in attribute `xyz_gen_heads`.

   - **`trees_to_lists()`**\
   Turns each `XYZGeneration` tree into lists representing Pauli paths, with the lists being appended to the attribute `xyz_pauli_paths`

   - **`xyz_tree_branching(cur_xyz_gen:XYZGeneration, pauli_paths_in_womb:List[List[PauliOperator]])`**\
   Recursively traverses all possible Pauli paths along an `XYZGeneration` tree and appends them to `xyz_pauli_paths`.

   - **`rn_to_z(first_op:PauliOperator)`**\
   A static method that replaces all "R"s and "N"s in the first `PauliOperator` of a path with "Z"s, in order to satisfy the second requirement to be a legal Pauli path.