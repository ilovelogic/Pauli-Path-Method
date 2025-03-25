from collections import defaultdict
from pauli_operator import PauliOperator
from layer import Layer
from pauli_path import PauliPath
from circuit import Circuit

def main():

    test_path = PauliPath(5, [3,2,3,2],[[(1,2),(3,4)],[(0,1),(2,3)],[(1,2),(3,4)],[(3,4)]])
    for key, sibs in test_path.all_ops[0].forward_sibs.items():
        print("Collection sorted based on {}".format(key))
        for op in sibs:
            print(op.operator)

    test_circuit = Circuit(6,9,17,[[(1,2),(3,4)],[(0,1),(2,3)],[(1,2),(3,4)],[(3,4)],[(1,2),(3,4)],[(0,1),(2,3)],[(1,2),(3,4)],[(3,4)]])
    for weight_combo in test_circuit.weight_combos:
        print(weight_combo)

if __name__ == "__main__":
    main()