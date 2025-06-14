from __future__ import annotations
import copy
from typing import List
from Pauli_Path_Method.Path_Generation.pauli_operator import PauliOperator

class XYZGeneration:
    """
    Initiates a list of all PauliOperator objects that map to the same next PauliOperator object list
    """
    def __init__(self, pauli_ops:List[PauliOperator],next_index:int, pauli_path:List[PauliOperator]):
        self.parent_ops = pauli_ops
        self.pauli_path = pauli_path

        if next_index == len(self.pauli_path):
            self.next_gen = None
        elif self.pauli_path[next_index].next_ops == None: # next_op is in the last layer of the Pauli path
            self.rp_to_z(self.pauli_path[next_index]) # must only use 'I's and 'Z's
        else:
            self.rnp_to_xyz(next_index)

    # We traverse the Pauli path from left to right. 
    # If we encounter an 'R' or 'N', we can do any of 'X', 'Y', and 'Z'.
    # If we encounter a 'P', we must use the Pauli at the same index in the prior layer 
    def rnp_to_xyz(self, next_index:int):
        r_pos_list = []
        n_pos_list = []

        self.fill_pos_lists(next_index, r_pos_list, n_pos_list)

        filled_n_list = []

        # We will fill in the 'N' positions with 'X', 'Y', and 'Z' in all possible combinations
        # resulting in 3 ** len(n_pos_list) configurations
        for i in range(3 ** len(n_pos_list)):
            filled_n_list.append(copy.deepcopy(self.pauli_path[next_index]))

        self.fill_in_pos(filled_n_list,n_pos_list,'X',0,0)
        self.fill_in_pos(filled_n_list,n_pos_list,'Y',0,3 ** (len(n_pos_list)-1))
        self.fill_in_pos(filled_n_list,n_pos_list,'Z',0,2 * (3 ** (len(n_pos_list)-1)))

        filled_rn_list = [[] for _ in range(len(filled_n_list))] # one grouping per n filling configuration

        for i in range(len(filled_n_list)):
            for j in range(3 ** len(r_pos_list)): # 3 options for every 'R' position
                filled_rn_list[i].append(copy.deepcopy(filled_n_list[i]))

            self.fill_in_pos(filled_rn_list[i],r_pos_list,'X',0,0)
            self.fill_in_pos(filled_rn_list[i],r_pos_list,'Y',0,3 ** (len(r_pos_list)-1))
            self.fill_in_pos(filled_rn_list[i],r_pos_list,'Z',0,2 * (3 ** (len(r_pos_list)-1)))
        
        self.next_gen = []
        for i in range(len(filled_n_list)):
            self.next_gen.append(XYZGeneration(filled_rn_list[i],next_index+1, self.pauli_path))

    def fill_pos_lists(self, next_index:int, r_pos_list: List[int], n_pos_list: List[int]):
        for i in range(len(self.pauli_path[next_index].operator)):
            if self.pauli_path[next_index].operator[i] == 'R':
                r_pos_list.append(i)
            elif self.pauli_path[next_index].operator[i] == 'N':
                carries = self.carries_to_the_end(next_index, i)
                if (carries == 1): # if the qubit remains a non-gate qubit to the last layer
                    self.pauli_path[next_index].operator[i] = 'Z' # to ensure that the last layer is all 'I's and 'Z's
                elif(carries == 2):
                    raise ValueError
                else:
                    n_pos_list.append(i)
            elif self.pauli_path[next_index].operator[i] == 'P': # the qubit is forced to be the same as the prior qubit
                self.pauli_path[next_index].operator[i] = self.parent_ops[0].operator[i]


    """
    This function recursively fills specified qubit positions of PauliOperators in
    filled_pos list with all possible combinations of 'X', 'Y', and 'Z'

    Args:
        self (XYZGeneration) :  The XYZGeneration instance, used to recursively call this function.
        filled_pos (List[PauliOperator]) : The list of PauliOperators to fill with all possible combinations
        of 'X', 'Y', and 'Z' at the positions in pos_list. This list is updated in place.
        pos_list (List[int]) : The list of qubit positions that we will recursively fill with
        non-identity Paulis.
        index (int) : The index of pos_list that contains the qubit position at which we want to 
        edit every PauliOperator in filled_pos in the range beginning with start.
        pauli (str) : The Pauli ('X', 'Y', or 'Z') to write to that given index .
        start (int) : The beginning of the range of PauliOperators in filled_pos that
        we should edit on this function call.

    Returns:
        void : Updates filled_pos to contain all PauliOperators whose qubits at the positions in pos_list
        contain all possible combinations of "X", "Y", and "Z".
    """
    def fill_in_pos(self,filled_pos:List[PauliOperator],pos_list:List[int], index:int, pauli:str, start:int):
        if index == len(pos_list): # We have filled all positions in pos_list
            return
        
        for i in range(start, start + (3 ** (len(pos_list)-index-1))):
            filled_pos[i].operator[pos_list[index]] = pauli
        
        self.fill_in_pos(filled_pos,pos_list,'X', index+1,start)
        self.fill_in_pos(filled_pos,pos_list,'Y', index+1,start+(3 ** (len(pos_list)-index-2)))
        self.fill_in_pos(filled_pos,pos_list,'Z', index+1,start+(2 * (3 ** (len(pos_list)-index-2))))

    # Checks if the non-gate qubit at index i in the Pauli operator at index pauli_path_index
    # in the Pauli path carries to the end of the Pauli path, i.e., if it remains a non-gate qubit
    # until the last layer of the Pauli path. In this case, it would be forced to 
    # always be a 'Z' to meet the all 'I's and 'Z' at last layer requirement
    def carries_to_the_end(self, pauli_path_index:int, i:int):
        if self.pauli_path[pauli_path_index].next_ops == None: # reached last layer
            return 1 # the non-gate qubit remain a non-gate qubit to the last layer
        elif (self.pauli_path[pauli_path_index].next_ops == []): # NEED TO CHECK IF THIS CAUSES ANY ISSUES
            return 2
        elif self.pauli_path[pauli_path_index].next_ops[0].operator[i] == 'P':
            return self.carries_to_the_end(pauli_path_index+1, i)
        return 0


    # The last Pauli operator in a Pauli path can only be a tensor of 'I's and 'Z's
    def rp_to_z(self, next_op:PauliOperator):
        for i in range(len(next_op.operator)):
            if next_op.operator[i] == 'R':
                next_op.operator[i] = 'Z'
            elif next_op.operator[i] == 'N':
                return 0 # there should be no 'N' in our last layer
            elif next_op.operator[i] == 'P':
                next_op.operator[i] = 'Z' # we set up our propagation to gurantee the prior of the last layer 
                # would have all Z's in non-gate qubit positions with Hamming weight
            self.next_gen = [XYZGeneration([next_op],len(self.pauli_path),self.pauli_path)] 
        return 1 # valid operator possible