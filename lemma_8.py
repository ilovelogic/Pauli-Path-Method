from typing import List


"""
For 2D grid architecture have a list of lists of lists
[
layer 1 (t=1): [first qubit position, second qubit position], etc.
layer 2 (t=2): [first qubit position, second qubit position], etc.
...
layer d (t=d): [first qubit position, second qubit position], etc.
]


Example:
[
[[0,1], [5,6], [3,7], [8,9], [10,11]],
[[0,4], [2,3], [5,9], [6,10], [7,11], [12,13], [14,15]],
[[1,2], [4,5], [6,7], [9,10], [8,12], [13,14]],
[[1,5], [2,6], [4,8], [9,13], [10,14], [11,15]]
]


Let q_{ij} be the qubit in the ith row and jth column of the circuit.
Let n be the total number of rows in the circuit and m be the total number of columns.
Assume we start indexing from 1.
Map each ij to m(i-1) + j-1.
This will result in unique identifying numbers from 0 to m(n-1)+m-1 = mn-1.
Valid qubit pairs are of the form [x,x+1] and [x, x+m].
The weight of layer l+1 can be at most the weight of layer l plus the number of gates between their layers 1


"""




def get_hamming_weights(d:int, l:int, n:int): # d = depth, l = upper bound on weight, n = number of qubits
   
    if (l > n*d):
        print("Invalid maximum Hamming weight")
        return []
   
    all_weight_combos = []
    enumerate_weights(all_weight_combos, [], l-d, d)
    return all_weight_combos


"""
This function determines all solutions to the equation w_1 + w_2 + ... + w_d = k, where w_i >= 1




It takes four arguments and recursively constructs our list of all weight combinations.




Args:
    weight_combos (List[List[int]]) : List of all possible weight combinations
    weight_list (List[int]) : List of the weight combination we're currently working with
    wiggle_room (int) : Number of layers that can have weight greater than 1
    num_layers_left (int) : Number of layers waiting for weight




Returns:
    void: The function appends to weight_combos, which affects the list in the calling function
"""


# Black box function
# w_1 + ... + w_d = k, w_i >= 1
# stars (d-k) and bars problem: look for already defined function
# tiny optimizations because of how layers map to each other


# ensure not bottleneck
def enumerate_weights(weight_combos:List[List[int]], weight_list:List[int], wiggle_room:int, num_layers_left:int):


   # base case: no more layers to add weight to
    if num_layers_left == 0:
        weight_combos.append(weight_list)
        return


    # base case: no more wiggle room
    if wiggle_room == 0:
        weight_list_copy = list(weight_list) # copy of weight_list to avoid overwriting
        for i in range(num_layers_left): # enumerates num_layers_left times, once for each remaining layer
            weight_list_copy.append(1) # all remaining layers get weight 1, since there's no wiggle room
        weight_combos.append(weight_list_copy)
        return


    # recursive case: still have wiggle room and layers to add weight to
    else:
        for i in range(wiggle_room+1):
            weight_list_copy = list(weight_list) # copy of weight_list to avoid overwriting
            weight_list_copy.append(i+1) # +1 since every weight is guaranteed to be least 1
            enumerate_weights(weight_combos, weight_list_copy, wiggle_room-i, num_layers_left-1)
        return








"""
This function determines all possible gate outputs/inputs of a new layer,
given the prior layer and the weights of the new layer and the prior layer.




It takes five arguments and uses helper function add_gate_input to build a list of all gate output/input combos




Args:
    prior_layer (str) : the sting representing the prior layer of the circuit,
    with each character representing a qubit
    prior_weight (int) : the Hamming weight of the prior layer
    num_gates (int) : the number of 2 qubit gates between the prior layer and the new layer
    this_weight (int) : the Hamming weight of the new layer
    odd_start (int) : 0 if the first 2 qubit gate takes as input the first two qubits,
    1 if it takes the second and third qubits
    n (int) : the number of qubits




Returns:
    List[str]: a list of all valid input/output combos for this layer, with each layer represented as a string
"""
def R_iterations(prior_layer:str, prior_weight:int, num_gates:int, this_weight:int, odd_start:int, n:int):


    if (num_gates > (n/2)):
        print("Invalid number of gates")
        return []
   
    pos_to_fill = [] # stores the positions in the prior_layer that need non-identity output
    for i in range(odd_start, num_gates*2, 2):
        if (prior_layer[i] == 'R' or prior_layer[i+1] == 'R'): # instead of i and i+1
        # enumerate over the index combinations ()
            pos_to_fill.append(i)


    num_RRs = this_weight - len(pos_to_fill) # number of RRs we can use to fill in the layer


    # distingush between rrs and irs
    if num_RRs < 0 or num_RRs > len(pos_to_fill):
        print("Could not satisfy user input with a valid layer")
        return [] # no way to make a valid layer, given the weights
    else:
        return add_gate_input([""], num_RRs, pos_to_fill, n)




def append_to_layers(layers:List[str], append_paulis:str):
    for i in range(len(layers)):
        layers[i] += append_paulis # only alters local copy
    return layers # returns all layers with apppend_paulis added


def add_gate_input(layers:List[str], num_RRs:int, pos_to_fill:List[int], n:int):


    if (len(pos_to_fill) == 0):
        while len(layers[0]) < n: # fill in the first layer with 'II' until we hit pos
            layers = append_to_layers(layers, 'I')
        return layers
   
    while len(layers[0]) < pos_to_fill[0]: # fill in the first layer with 'II' until we hit pos
        layers = append_to_layers(layers, 'I') # maybe add 1 I if odd start and then add in chunks of II after that to optimize?
   
    pos_to_fill.pop(0) # remove the position we just handled the 'II's up to


    if (len(pos_to_fill)+1 == num_RRs): # no more wiggle room, we must fill all remaining gate inputs with RR    
        rr_layers = append_to_layers(list(layers), 'RR')
        return add_gate_input(rr_layers, num_RRs-1, list(pos_to_fill), n) # one less RR to use,
        # next call will handle adding the next RR
    else:
        ir_layers = append_to_layers(list(layers), 'IR') # copy of layers with 'IR' added to all strs


        ri_layers = append_to_layers(list(layers), 'RI') # copy of layers with 'RI' added to all strs


        if (num_RRs == 0):
            return add_gate_input(ir_layers, num_RRs, list(pos_to_fill), n) + add_gate_input(ri_layers, num_RRs, list(pos_to_fill), n)


        rr_layers = append_to_layers(list(layers), 'RR') # copy of layers with 'RR' added to all strs


        return add_gate_input(ir_layers, num_RRs, list(pos_to_fill), n) + add_gate_input(ri_layers, num_RRs, list(pos_to_fill), n) + add_gate_input(rr_layers, num_RRs-1, list(pos_to_fill), n)
        #return ir_layers + ri_layers + rr_layers # append all possible combinations together


def main():
    # testing weight enumeration
    all_weight_combos = get_hamming_weights(3, 6, 4)
    print()
    print("Weight combinations for d=3, l=6, n=4:")
    for sublist in all_weight_combos:
        for element in sublist:
            print(element, end=", ")
        print() # move to the next line after each sublist
    print()


    # testing layer propagation
    print("R_iterations for prior layer IRIRRRRI, new weight = 6")
    #R_iterations(prior_layer:str, prior_weight:int, num_gates:int, this_weight:int, odd_start:int, n:int)
    all_layer_props = R_iterations("IRIRRRRI", 5, 4, 6, 0, 8)
    for element in all_layer_props:
        print(element)
       


if __name__ == "__main__":
    main()