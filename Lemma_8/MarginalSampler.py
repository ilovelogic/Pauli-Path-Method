from collections import defaultdict
import numpy as np
from Pauli_Amplitude.edited_pauli_amp import compute_marginal_noisy_fourier


class MarginalSampler:
    def __init__(self, C, sib_op_heads, n_qubits, gamma, num_samples=16):
        self.C = C
        self.sib_op_heads = sib_op_heads
        self.n = n_qubits
        self.gamma = gamma
        self.num_samples = num_samples

        # Sample immediately and store normalized distribution
        #self.sampled_probs = self.sample_many(1)  # was 2000

        self.sampled_probs = self.sample_many(self.num_samples)

    def sample(self):
        return self.sample_marginal_bits()

    def sample_many(self, num_samples):
        counts = defaultdict(int)
        for _ in range(num_samples):
            bitstring = self.sample()
            counts[bitstring] += 1

        # Normalize to get probability distribution
        total = sum(counts.values())
        return {x: c / total for x, c in counts.items()}

    def sample_marginal_bits(self):
        """
        Samples a full bitstring x ∈ {0,1}^n using bit-by-bit marginal sampling.
        """
        fixed_bits = {}

        for i in range(self.n):
            fb0 = fixed_bits.copy()
            fb0[i] = '0'
            fb1 = fixed_bits.copy()
            fb1[i] = '1'

            p0 = compute_marginal_noisy_fourier(self.C, self.sib_op_heads, fb0, self.n, self.gamma)
            p1 = compute_marginal_noisy_fourier(self.C, self.sib_op_heads, fb1, self.n, self.gamma)

            total = p0 + p1
            if total == 0:
                prob0 = 0.5
            else:
                prob0 = p0 / total

            bit = '0' if np.random.rand() < prob0 else '1'

            #reversed_fb0 = {self.n - 1 - i: b for i, b in fb0.items()}
            #p0 = compute_marginal_noisy_fourier(self.C, self.sib_op_heads,fb0, self.n, self.gamma)

            #prob0 = min(max(p0, 0.0), 1.0)  # ensure in [0, 1] range
            #print('prob0', prob0)
            #print("bit:", bit)
            #print(f"[DEBUG] Qubit {i}: p(0|...) = {p0:.4e}, p(1|...) = {p1:.4e}, chosen bit = {bit}")
            fixed_bits[i] = bit
        #print('fixed_bits', fixed_bits)
        return ''.join(fixed_bits[i] for i in reversed(range(self.n)))
