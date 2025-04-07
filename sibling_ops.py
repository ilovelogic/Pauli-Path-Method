from __future__ import annotations
import copy
from typing import List
from pauli_operator import PauliOperator
from pauli_path_trav import PauliPathTrav

class SiblingOps:
    """
    Initiates a list of all PauliOperator objects that map to the same next PauliOperator object list
    """
    def __init__(self, pauli_ops:List[PauliOperator], next_op:PauliOperator):
        
        self.pauli_ops = pauli_ops
        self.rnp_to_xyz(next_op)

     # Must traverse each trav from left to right. 
    # If we encounter an 'R' or 'N', we can do any of 'X', 'Y', and 'Z'.
    # If we encounter a 'P', we must use the Pauli at the same index in the prior layer 
    def rnp_to_xyz(self, next_op:PauliOperator):
        r_pos_list = []
        n_pos_list = []
        filled_n_list = []
        for i in range(len(next_op.operator)):
            if next_op.operator[i] == 'R':
                r_pos_list.append(i)
            elif next_op.operator[i] == 'N':
                if (self.carries_to_the_end(next_op.operator[0],i)): # if the qubit remains a non-gate qubit to the last layer
                    next_op.operator[i] = 'Z' # to ensure that the last layer is all 'I's and 'Z's
                else:
                    n_pos_list.append(i)
            elif next_op.operator[i] == 'P': # the qubit is forced to be the same as the prior qubit
                next_op.operator[i] = self.pauli_ops[0].operator[i]

        for i in range(3 ** len(n_pos_list)):
            filled_n_list.append(copy.deepcopy(next_op))
            
        self.fill_in_pos(filled_n_list,n_pos_list,0,0)

        next_ops = []*len(filled_n_list)

        for i in range(len(filled_n_list)):
            for j in range(3 ** len(r_pos_list)):
                next_ops[i].append(copy.deepcopy(filled_n_list[i]))
            self.fill_in_pos(filled_n_list[i],r_pos_list,'X',0,0)
            self.fill_in_pos(filled_n_list[i],r_pos_list,'Y',0,3 ** (len(r_pos_list)-1))
            self.fill_in_pos(filled_n_list[i],r_pos_list,'Z',0,2 * (3 ** (len(r_pos_list)-1)))

    def fill_in_pos(self,filled_pos:List[PauliOperator],pos_list:List[int], pauli:str, index:int, start:int):
        if index == len(pos_list):
            return
        for i in range(start, start + (3 ** (len(pos_list)-index-1))):
            filled_pos[i][pos_list[index]] = pauli
        
        self.fill_in_pos(filled_pos,pos_list,'X', index+1,start)
        self.fill_in_pos(filled_pos,pos_list,'Y', index+1,start+(3 ** (len(pos_list)-index-2)))
        self.fill_in_pos(filled_pos,pos_list,'Z', index+1,start+(2 * (3 ** (len(pos_list)-index-2))))

    def carries_to_the_end(self, cur_op:PauliOperator, i:int):
        if cur_op.next_ops == None: # reached last layer
            return True # the non-gate qubit remain a non-gate qubit to the last layer
        elif cur_op.next_ops[0].operator[i] == 'P' or cur_op.next_ops[0].operator[i] == 'N':
            return self.carries_to_the_end(cur_op.next_ops[0], i)
        return False


    # The first and last Pauli operator in a Pauli path can only be a tensor of 'I's and 'Z's
    def rn_to_z(self, cur_op:PauliOperator, prior_op:PauliOperator=None):
        for i in range(len(cur_op.operator)):
            if cur_op.operator[i] == 'R' or cur_op.operator[i] == 'N':
                cur_op.operator[i] = 'Z'
            if self.operator[i] == 'P' and prior_op != None: # only possible when we are at the last (not first) op
                if prior_op[i] == 'Z':
                    cur_op.operator[i] = 'Z'
                else:
                    return 0 # indicates the need to remove the path
        return 1 # valid path possible