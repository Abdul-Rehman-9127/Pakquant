"""PakQuant - quantum solver: QUBO -> Ising Hamiltonian -> QAOA (Qiskit).

Pipeline
  1. QUBO   : minimise  -sum(score_i * x_i) + PENALTY * (sum(x_i) - 3)^2
  2. Ising  : substitute x_i = (1 - Z_i)/2  ->  SparsePauliOp of I / Z / ZZ terms
  3. QAOA   : qaoa_ansatz(H, reps=p), angles tuned by COBYLA on <H>
  4. Sample : 1024 shots from the tuned circuit

Bit convention (identical to brute.py): qubit i  <->  institution i
(0 = ICU ... 4 = mall).  Qiskit PRINTS bitstrings as q4 q3 q2 q1 q0, so the
optimal answer ICU+dialysis+water is printed by Qiskit as "00111" - the same
string brute.py prints.  We also show it in institution order (x1..x5).

Technical notes
  * Scores are divided by 100 so the Hamiltonian has O(1) coefficients.  QAOA
    angles are periodic; with raw coefficients (e.g. a 500 penalty) the gamma
    landscape oscillates wildly and COBYLA gets stuck.
  * PENALTY = 1.0 (= 100 lives, the largest single score) is a safe bound:
    any A >= max(w) guarantees that no constraint-violating state can beat the
    best valid one.  (For this instance the true minimum is only A > 0.05, so
    1.0 is conservative.)
  * COBYLA is local, so we use several seeded random restarts and keep the
    run with the LOWEST ENERGY (we never look at the answer to pick a run).
"""
import numpy as np
from scipy.optimize import minimize
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.circuit.library import qaoa_ansatz
from qiskit.primitives import StatevectorEstimator, StatevectorSampler

# ---------------- problem ----------------
names = ["ICU", "dialysis", "water", "school", "mall"]
scores = np.array([100, 80, 40, 5, 1])
n, K = 5, 3
SCALE = 100.0
w = scores / SCALE            # normalised scores
PENALTY = 1.0                 # in normalised units (= 100 lives)

P_LAYERS = 2                  # headline circuit depth
SWEEP = [1, 2, 3]             # depths shown in the comparison table
RESTARTS, MAXITER = 15, 200
SHOTS, SEED = 1024, 7


def cost(bits):
    """QUBO objective (normalised units) for a list of n bits."""
    return -sum(w[i] * bits[i] for i in range(n)) + PENALTY * (sum(bits) - K) ** 2


def build_hamiltonian():
    """Ising Hamiltonian from the QUBO via x_i = (1 - Z_i)/2.
    Qiskit Pauli labels are little-endian: qubit i is label[n-1-i]."""
    terms = []
    const = -sum(w) / 2 + PENALTY * (n / 2 - K) ** 2 + PENALTY * n / 4
    for i in range(n):
        lab = ["I"] * n
        lab[n - 1 - i] = "Z"
        terms.append(("".join(lab), w[i] / 2 - PENALTY * (n / 2 - K)))
    for i in range(n):
        for j in range(i + 1, n):
            lab = ["I"] * n
            lab[n - 1 - i] = lab[n - 1 - j] = "Z"
            terms.append(("".join(lab), PENALTY / 2))
    return SparsePauliOp.from_list(terms + [("I" * n, const)])


def decode(qiskit_bits):
    """'00111' (Qiskit order q4..q0) -> (x1..x5 string, names, lives score, valid?)"""
    x = [int(qiskit_bits[n - 1 - i]) for i in range(n)]
    chosen = [names[i] for i in range(n) if x[i]]
    return "".join(map(str, x)), chosen, int(sum(scores[i] * x[i] for i in range(n))), sum(x) == K


def run_qaoa(H, p, seed):
    ansatz = qaoa_ansatz(H, reps=p)
    est = StatevectorEstimator()
    energy = lambda th: float(est.run([(ansatz, H, th)]).result()[0].data.evs)
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(RESTARTS):
        r = minimize(energy, rng.uniform(0, np.pi, 2 * p), method="COBYLA",
                     options={"maxiter": MAXITER})
        if best is None or r.fun < best.fun:
            best = r
    return ansatz, best


if __name__ == "__main__":
    H = build_hamiltonian()
    print("Ising Hamiltonian:\n", H, "\n")

    # --- sanity check: Hamiltonian diagonal must equal the QUBO cost ---
    diag = np.real(np.diag(H.to_matrix()))
    ok = all(abs(diag[k] - cost([(k >> i) & 1 for i in range(n)])) < 1e-9 for k in range(2 ** n))
    print(f"Sanity check, all 32 bitstrings: Hamiltonian == QUBO cost -> {ok}")
    assert ok

    # --- depth sweep (exact output probabilities, no sampling noise) ---
    valid_idx = [k for k in range(2 ** n) if bin(k).count("1") == K]
    best_k = max(valid_idx, key=lambda k: sum(scores[i] for i in range(n) if (k >> i) & 1))
    opt_score = int(sum(scores[i] for i in range(n) if (best_k >> i) & 1))  # exact optimum, as in brute.py
    print("\nDepth sweep (best of %d COBYLA restarts, chosen by lowest energy):" % RESTARTS)
    print(f"  {'p':>2} {'energy':>9} {'P(optimal)':>11} {'P(valid)':>9}  top state (Qiskit order)")
    results = {}
    for p in SWEEP:
        ansatz, res = run_qaoa(H, p, SEED)
        probs = Statevector(ansatz.assign_parameters(res.x)).probabilities()
        results[p] = (ansatz, res)
        print(f"  {p:>2} {res.fun:>9.4f} {probs[best_k]:>11.1%} {probs[valid_idx].sum():>9.1%}  "
              f"{int(np.argmax(probs)):05b}")

    # --- headline run: p = 2, sampled 1024 shots ---
    ansatz, res = results[P_LAYERS]
    print(f"\nHeadline run: p={P_LAYERS}, final energy <H> = {res.fun:.4f}")
    print("Best angles  beta =", np.round(res.x[:P_LAYERS], 3), " gamma =", np.round(res.x[P_LAYERS:], 3),
          "(qaoa_ansatz orders its parameters beta first, then gamma)")
    circ = ansatz.assign_parameters(res.x)
    circ.measure_all()
    counts = StatevectorSampler(seed=SEED).run([circ], shots=SHOTS).result()[0].data.meas.get_counts()

    print(f"\nTop 10 of {SHOTS} shots   (Qiskit string = q4..q0 | x1..x5 = ICU..mall)")
    for bits, c in sorted(counts.items(), key=lambda kv: -kv[1])[:10]:
        x, chosen, sc, valid = decode(bits)
        tag = f"score {sc}" if valid else "INVALID (not exactly 3)"
        print(f"  {bits} | {x}  {c:4d} ({c / SHOTS:5.1%})  {' + '.join(chosen):<32} {tag}")

    top_bits, top_c = max(counts.items(), key=lambda kv: kv[1])
    x, chosen, sc, valid = decode(top_bits)
    n_valid = sum(c for b, c in counts.items() if decode(b)[3])
    print(f"\nMost frequent outcome: {' + '.join(chosen)} (score {sc}) "
          f"-> {'MATCHES' if sc == opt_score else 'does NOT match'} the brute-force optimum ({opt_score})")
    print(f"Optimal-answer frequency: {top_c}/{SHOTS} = {top_c / SHOTS:.1%}  "
          f"(random guessing over 32 states: 3.1%, uniform over the 10 valid: 10%)")
    print(f"Shots that were valid (exactly 3): {n_valid}/{SHOTS} = {n_valid / SHOTS:.1%}")
