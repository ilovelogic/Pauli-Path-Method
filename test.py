import time
from collections import defaultdict
from pauli_operator import PauliOperator
from layer import Layer
from pauli_path import PauliPath
from circuit import Circuit

def main():

    print(time.time())
    test_path = PauliPath(5, [3,2,3,2],[[(1,2),(3,4)],[(0,1),(2,3)],[(1,2),(3,4)],[(3,4)]])
    for key, sibs in test_path.all_ops[0].forward_sibs.items():
        print("Collection sorted based on {}".format(key))
        for op in sibs:
            print(op.operator)
    #test_circuit = Circuit(5,4,8,[[(1,2),(3,4)],[(0,1),(2,3)],[(1,2),(3,4)]])
    test_gate_pos = [[]]*6
    for i in range(3):
        for j in range(0,12,2):
            test_gate_pos[i].append((j,j+1))
        for j in range(1,13,2):
            test_gate_pos[i].append((j,j+1))
    test_circuit = Circuit(15,7,21,test_gate_pos)
    for pauli_path in test_circuit.pauli_paths:
        for sibs in pauli_path.layers[5].forward_sibs.values():
            for op in sibs:
                print(op.operator)
    pauli = PauliOperator(["I","R","R","I","R","R"])
    pauli.r_to_xyz()
    for pauli_op in pauli.xyz_paulis:
        print(pauli_op)
    print(time.time())

if __name__ == "__main__":
    main()