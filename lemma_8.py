from typing import List

def list_allocs(num_p:int, num_w:int):
    list_alloc = [[0 for _ in range(num_w+1)] for _ in range(num_p+1)]

    # base cases
    for w in range(0, num_w):
        list_alloc[0][w] = 0 # no gates, no way to have a layer

    for p in range(0, num_p): # no weight, only one way to have a layer
        list_alloc[p][0] = 1
    
    for i in range(1, num_p+1):
        for j in range(1,num_w+1):   #     IR                     RI                     RR
            if (j > i):
                list_alloc[i][j] = list_alloc[i-1][j-1] + list_alloc[i-1][j-1] + list_alloc[i-1][j-2]
            else:
                list_alloc[i][j] = list_alloc[i-1][j-1] + list_alloc[i-1][j-1]
            print(f"list_alloc[{i}][{j}] = {list_alloc[i][j]}")
        
    return list_alloc


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
For 2D grid architecture have a list of lists of tuples
[
layer 1 (t=1): [first qubit position, second qubit position], etc.
layer 2 (t=2): [first qubit position, second qubit position], etc.
...
layer d (t=d): [first qubit position, second qubit position], etc.
]

Example:
[
[(0,1), (5,6), (3,7), (8,9), (10,11)],
[(0,4), (2,3), (5,9), (6,10), (7,11), (12,13), (14,15)],
[(1,2), (4,5), (6,7), (9,10), (8,12), (13,14)],
[(1,5), (2,6), (4,8), (9,13), (10,14), (11,15)]
]

Let q_{ij} be the qubit in the ith row and jth column of the circuit. 
Let n be the total number of rows in the circuit and m be the total number of columns.
Assume we start indexing from 1.
Map each ij to m(i-1) + j-1. 
This will result in unique identifying numbers from 0 to m(n-1)+m-1 = mn-1.
Valid qubit pairs are of the form [x,x+1] and [x, x+m].
The weight of layer l+1 can be at most the weight of layer l plus the number of gates between their layers

"""

layers = [] # global variable

"""
This function determines all possible gate outputs/inputs of a new layer,
given the prior layer and the weights of the new layer and the prior layer.


It takes five arguments and uses helper function add_gate_input to build a list of all gate output/input combos


Args:
    prior_layer (str) : the sting representing the prior layer of the circuit,
    with each character representing a qubit
    this_weight (int) : the Hamming weight of the new layer
    n (int) : the number of qubits
    pos_list (List[tuple]) : each tuple contains the indices of the two inputs to one of the gates


Returns:
    List[str]: a list of all valid input/output combos for this layer, with each layer represented as a string
"""
def R_iterations(prior_layer:str, this_weight:int, n:int, pos_list:List[tuple]):

    if (len(pos_list) > (n/2)):
        print("Invalid number of gates")
        return

    layer_str = []

    for i in range(n):
        layer_str.append('I')
    
    pos_to_fill = [] # stores the positions in the prior_layer that need non-identity output
    for ind1, ind2 in pos_list:
        if (prior_layer[ind1] == 'R' or prior_layer[ind2] == 'R'): # instead of i and i+1
            # enumerate over the index combinations
            pos_to_fill.append((ind1,ind2))

    for i in range(3**(len(pos_to_fill))): # 3^{len(pos_to_fill)} possible selections of IR, RI, and RR
        layers.append(list(layer_str)) # copies layer_str

    num_RRs = this_weight - len(pos_to_fill) # number of RRs we can use to fill in the layer

    # distingush between rrs and irs
    if num_RRs < 0 or num_RRs > len(pos_to_fill):
        print("Could not satisfy user input with a valid layer")
        # no way to make a valid layer, given the weights
    else:
        add_gate_input(num_RRs, pos_to_fill, 0, 3**(len(pos_to_fill)))

def append_to_layers(indices:tuple, strs:tuple, r_start:int, r_end:int):
    ind1, ind2 = indices
    str1, str2 = strs
    for i in range(r_start, r_end):
        layers[i][ind1] = str1
        layers[i][ind2] = str2

# make layers a global variable

def add_gate_input(num_RRs:int, pos_to_fill:List[tuple], r_start:int, r_end:int):

    if (len(pos_to_fill) == 0):
        return
    
    if (len(pos_to_fill) == num_RRs): # no more wiggle room, we must fill all remaining gate inputs with RR     
        cur_pos = pos_to_fill.pop(0)
        rr_start = r_start+2*(3**len(pos_to_fill))
        rr_end = r_start+3*(3**len(pos_to_fill))
        append_to_layers(cur_pos, ('R','R'), rr_start, rr_end) # copy of layers with 'RR' added to all strs
        add_gate_input(num_RRs-1, list(pos_to_fill), rr_start, rr_end) # one less RR to use,
        # next call will handle adding the next RR
        return

    else:
        # IR -> r_start += 0*3^{len(pos_to_fill)-1}
        # RI -> r_start +=  1*3^{len(pos_to_fill)-1}
        # RR -> r_start += 2*3^{len(pos_to_fill)-1}
        #       r_end = r_start + 1*3^{len(pos_to_fill)-1} // not inclusive
        # base 3 string, length len(pos_to_fill)?
        cur_pos = pos_to_fill.pop(0)

        ir_start = r_start
        ir_end = r_start+3**len(pos_to_fill)
        append_to_layers(cur_pos, ('I','R'), ir_start, ir_end) # copy of layers with 'IR' added to all strs

        ri_start = r_start+3**len(pos_to_fill)
        ri_end = r_start+2*(3**len(pos_to_fill))
        append_to_layers(cur_pos, ('R','I'), ri_start, ri_end) # copy of layers with 'RI' added to all strs

        if (num_RRs != 0):
            rr_start = r_start+2*(3**len(pos_to_fill))
            rr_end = r_start+3*(3**len(pos_to_fill))
            append_to_layers(cur_pos, ('R','R'), rr_start, rr_end) # copy of layers with 'RR' added to all strs
            add_gate_input(num_RRs-1, list(pos_to_fill), rr_start, rr_end) 

        add_gate_input(num_RRs, list(pos_to_fill), ir_start, ir_end) 
        add_gate_input(num_RRs, list(pos_to_fill), ri_start, ri_end) 
        return

def main():
    list_alloc = list_allocs(7,8)
    # testing weight enumeration
    #all_weight_combos = get_hamming_weights(3, 6, 4)
    #print()
    #print("Weight combinations for d=3, l=6, n=4:")
    #for sublist in all_weight_combos:
        #for element in sublist:
            #print(element, end=", ")
        #print() # move to the next line after each sublist
    #print()

    # testing layer propagation
    #print("R_iterations for prior layer IRIIR, new weight = 3")
    #R_iterations(prior_layer:str, this_weight:int, n:int)
    #R_iterations("IRIRI", 3, 5, [(0,1),(2,3)])
    #for lil_list in layers:
        #str = ''
        #for c in lil_list:
            #str += c
        #print(str)

if __name__ == "__main__":
    main()