from collections import defaultdict
from pauli_operator import PauliOperator
from layer import LayerOperators
from pauli_path import PauliPath
from circuit import Circuit

def main():

    # 3 qubit TEST CASE!!!

    
    test_path = PauliPath(5, [3,2,3],[[(1,2),(3,4)],[(0,1),(2,3)],[(1,2),(3,4)]])
    for key, sibs in test_path.all_ops[1].backward_sibs_list.items():
        print("Collection sorted based on {}".format(key))
        for op in sibs:
            print(op.operator)

if __name__ == "__main__":
    main()