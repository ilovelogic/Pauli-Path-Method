from __future__ import annotations
import copy
from typing import List
from pauli_operator import PauliOperator

class SiblingOps:
    """
    Initiates a list of all PauliOperator objects that map to the same next PauliOperator object list
    """
    def __init__(self, pauli_ops:List[PauliOperator],next_index:int, pauli_path:List[PauliOperator]):
        self.pauli_ops = pauli_ops
        self.pauli_path = pauli_path

        if next_index == len(self.pauli_path):
            self.next_sibs = None
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

        for i in range(3 ** len(n_pos_list)):
            filled_n_list.append(copy.deepcopy(self.pauli_path[next_index]))

        self.fill_in_pos(filled_n_list,n_pos_list,'X',0,0)
        self.fill_in_pos(filled_n_list,n_pos_list,'Y',0,3 ** (len(n_pos_list)-1))
        self.fill_in_pos(filled_n_list,n_pos_list,'Z',0,2 * (3 ** (len(n_pos_list)-1)))

        filled_rn_list = [[] for _ in range(len(filled_n_list))] # one grouping per n filling configuration

        for i in range(len(filled_n_list)):
            for j in range(3 ** len(r_pos_list)):
                filled_rn_list[i].append(copy.deepcopy(filled_n_list[i]))

            self.fill_in_pos(filled_rn_list[i],r_pos_list,'X',0,0)
            self.fill_in_pos(filled_rn_list[i],r_pos_list,'Y',0,3 ** (len(r_pos_list)-1))
            self.fill_in_pos(filled_rn_list[i],r_pos_list,'Z',0,2 * (3 ** (len(r_pos_list)-1)))
        
        self.next_sibs = []
        for i in range(len(filled_n_list)):
            self.next_sibs.append(SiblingOps(filled_rn_list[i],next_index+1, self.pauli_path))

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
                self.pauli_path[next_index].operator[i] = self.pauli_ops[0].operator[i]


    def fill_in_pos(self,filled_pos:List[PauliOperator],pos_list:List[int], pauli:str, index:int, start:int):
        if index == len(pos_list): # index is the position of the index in the pos_list
            # at which we want to edit the pauli in filled_pos's operators in the range beginning with start
            return
        for i in range(start, start + (3 ** (len(pos_list)-index-1))):
            filled_pos[i].operator[pos_list[index]] = pauli
        
        self.fill_in_pos(filled_pos,pos_list,'X', index+1,start)
        self.fill_in_pos(filled_pos,pos_list,'Y', index+1,start+(3 ** (len(pos_list)-index-2)))
        self.fill_in_pos(filled_pos,pos_list,'Z', index+1,start+(2 * (3 ** (len(pos_list)-index-2))))

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
            self.next_sibs = [SiblingOps([next_op],len(self.pauli_path),self.pauli_path)] 
        return 1 # valid operator possible