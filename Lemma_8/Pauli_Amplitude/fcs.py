

class FCalculator:
    def __init__(self, C, s, x):
        self.C = C          # Circuit layers
        self.s = s          # Pauli path (list of strings)
        self.x = x          # Output state (binary string)
        self.num_qubits = len(s[0]) if s else 0


    def calculate_f(self):
        # Check valid Pauli path length
        if len(self.s) != len(self.C) + 1:
            return 0.0
        
        # Check final Pauli is I/Z only
        final_pauli = self.s[-1]
        if any(p not in {'I', 'Z'} for p in final_pauli):
            return 0.0
        
        # Compute measurement factor (-1)^(sum x_i * z_i)
        measurement_factor = 1.0
        for i, (pauli, bit) in enumerate(zip(final_pauli, self.x)):
            if pauli == 'Z' and bit == '1':
                measurement_factor *= -1
        
        # Track total phase from gate transitions
        total_phase = 1.0
        for layer_idx, layer in enumerate(self.C):
            current_pauli = self.s[layer_idx]
            next_pauli = self.s[layer_idx + 1]
            
            for gate, qubits in layer:
                # Extract input/output Paulis for the gate's qubits
                input_sub = ''.join([current_pauli[q] for q in qubits])
                output_sub = ''.join([next_pauli[q] for q in qubits])
                
                # Compute phase contribution for this gate
                phase = self._get_gate_phase(gate, input_sub, output_sub)
                if phase is None:
                    return 0.0  # Invalid transition
                total_phase *= phase
        
        return total_phase * measurement_factor

    def _get_gate_phase(self, gate, input_pauli, output_pauli):
        # Placeholder: Implement phase logic for specific gates
        # Example for CNOT (assumes gate is CNOT and qubits = [control, target])
        if "CNOT" in str(gate):
            control, target = 0, 1  # Assuming qubits are ordered [control, target]
            c_in, t_in = input_pauli[control], input_pauli[target]
            c_out, t_out = output_pauli[control], output_pauli[target]
            
            # CNOT transition rules
            if c_in == 'X':
                expected_t_out = 'X' if t_in == 'I' else 'I' if t_in == 'X' else None
            elif c_in == 'Z':
                expected_t_out = t_in
            else:
                # Handle Y or other cases (simplified)
                return None  # Extend for full support
            
            if (t_out != expected_t_out) or (c_out != c_in):
                return None  # Invalid transition
            
            # Phase calculation (simplified; extend for full accuracy)
            return 1.0 if c_in != 'Y' else -1.0
        
        # Extend for other gates (e.g., Hadamard, Phase)
        return None  # Unsupported gate