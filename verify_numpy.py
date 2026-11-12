"""PakQuant - independent cross-check of q1.py WITHOUT Qiskit.

Re-implements the same QAOA from scratch with NumPy (uniform superposition,
then p x [ phase exp(-i*gamma*C) , mixer RX(2*beta) on every qubit ]), using the
same QUBO, the same seeded COBYLA restarts and the same parameter order that
qaoa_ansatz uses (all betas first, then all gammas).

If the depth-sweep table below agrees with q1.py's output, the Qiskit pipeline is
doing what the maths says.  Expect p=1 and p=2 to match q1.py exactly; the p=3 row
can differ in the last digit (tiny floating-point differences change COBYLA's path).

Run:  python verify_numpy.py      (needs only numpy and scipy)
"""
import numpy as np
from scipy.optimize import minimize

names = ["ICU", "dialysis", "water", "school", "mall"]
scores = np.array([100, 80, 40, 5, 1])
n, K = 5, 3
w, PENALTY = scores / 100.0, 1.0
RESTARTS, MAXITER, SEED = 15, 200, 7


def cost(k):  # same QUBO as q1.py; bit i of k <-> institution i
    b = [(k >> i) & 1 for i in range(n)]
    return -sum(w[i] * b[i] for i in range(n)) + PENALTY * (sum(b) - K) ** 2


diag = np.array([cost(k) for k in range(2 ** n)])
valid = [k for k in range(2 ** n) if bin(k).count("1") == K]
best_k = max(valid, key=lambda k: sum(scores[i] for i in range(n) if (k >> i) & 1))


def mixer(psi, beta):  # RX(2*beta) on every qubit
    c, s = np.cos(beta), -1j * np.sin(beta)
    psi = psi.reshape([2] * n)
    for ax in range(n):
        a, b = np.take(psi, 0, axis=ax), np.take(psi, 1, axis=ax)
        psi = np.stack([c * a + s * b, s * a + c * b], axis=ax)
    return psi.reshape(-1)


def probs(theta, p):
    betas, gammas = theta[:p], theta[p:]
    psi = np.ones(2 ** n, complex) / np.sqrt(2 ** n)
    for l in range(p):
        psi = mixer(psi * np.exp(-1j * gammas[l] * diag), betas[l])
    return np.abs(psi) ** 2


print(f"Exact optimum (enumeration): {' + '.join(names[i] for i in range(n) if (best_k >> i) & 1)}, "
      f"Qiskit-order bits {best_k:05b}\n")
print(f"{'p':>2} {'energy':>9} {'P(optimal)':>11} {'P(valid)':>9}  top state")
for p in (1, 2, 3):
    rng, best = np.random.default_rng(SEED), None
    for _ in range(RESTARTS):
        r = minimize(lambda t: float(probs(t, p) @ diag), rng.uniform(0, np.pi, 2 * p),
                     method="COBYLA", options={"maxiter": MAXITER})
        if best is None or r.fun < best.fun:
            best = r
    pr = probs(best.x, p)
    print(f"{p:>2} {best.fun:>9.4f} {pr[best_k]:>11.1%} {pr[valid].sum():>9.1%}  {int(np.argmax(pr)):05b}")
