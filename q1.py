from qiskit.quantum_info import SparsePauliOp
from qiskit import QuantumCircuit
from qiskit.circuit.library import qaoa_ansatz
from qiskit.primitives import StatevectorEstimator, StatevectorSampler
from scipy.optimize import minimize
import numpy as np


scores = [100, 80, 40, 5, 1]
names = ["ICU", "dialysis", "water", "school", "mall"]
n = 5
k_target = 3
penalty = 500
def cost(bits):
    chosen_score = 0
    for i in range(n):
        chosen_score = chosen_score + scores[i] * bits[i]
    count = sum(bits)
    penalty_term = penalty * (count - k_target) ** 2
    return -chosen_score + penalty_term



pauli_terms = []


constant = 0
for i in range(n):
    constant += -scores[i] / 2
constant += penalty * (n/2 - k_target)**2 + penalty * n/4


for i in range(n):
   
    coef = scores[i] / 2
   
    coef += -penalty * (n/2 - k_target)

    label = ["I"] * n
    label[n-1-i] = "Z"
    pauli_terms.append(("".join(label), coef))


for i in range(n):
    for j in range(i+1, n):
        label = ["I"] * n
        label[n-1-i] = "Z"
        label[n-1-j] = "Z"
        pauli_terms.append(("".join(label), penalty / 2))


H = SparsePauliOp.from_list(pauli_terms + [("I"*n, constant)])
print(H)
import numpy as np
H_matrix = H.to_matrix()
diag = np.real(np.diag(H_matrix))

print("Checking a few bitstrings:")
for bits in [[1,1,1,0,0], [1,1,1,1,1], [1,1,0,0,0], [1,0,1,0,1]]:
    k_qiskit = 0
    for i in range(n):
        if bits[i] == 1:
            k_qiskit += 2**i
    print(f"  bits={bits}  cost()={cost(bits)}  H diagonal={diag[k_qiskit]}")
   
from scipy.optimize import minimize

p = 2
ansatz = qaoa_ansatz(H, reps=p)
estimator = StatevectorEstimator()

def expectation(params):
    job = estimator.run([(ansatz, H, params)])
    return job.result()[0].data.evs

initial_params = np.array([0.5] * (2 * p))
print("Starting QAOA optimization...")
result = minimize(expectation, initial_params, method='COBYLA', options={'maxiter': 100})
print("Done. Final energy:", result.fun)
print("Best angles:", result.x)
tuned_circuit = ansatz.assign_parameters(result.x)
tuned_circuit.measure_all()

sampler = StatevectorSampler()
samples = sampler.run([tuned_circuit], shots=1024).result()
counts = samples[0].data.meas.get_counts()

print("\nTop 10 measurement results:")
sorted_counts = sorted(counts.items(), key=lambda x: -x[1])
for bitstring, count in sorted_counts[:10]:
    print(f"  {bitstring}: {count} times")