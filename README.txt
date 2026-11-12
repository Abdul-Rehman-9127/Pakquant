PakQuant — Qiskit Fall Fest 2026 Hackathon
============================================

Team: Muhammad Abdul Rehman (26i-0526)
      Muhammad Shehryar Malik
      Muhammad Taha Ikram
Institution: FAST-NUCES Islamabad

Problem
-------
A grid operator during a power shortage can only keep 3 of 5 critical
institutions powered this hour. Each institution is scored by lives at
risk if left without power: ICU=100, dialysis center=80, water pumping
station=40, school=5, mall=1. The goal is to choose the 3 institutions
that minimize harm.

Files
-----
brute.py  - classical exact baseline (brute-force, enumerates all 32
            combinations of 5 institutions, keeps the best valid one)
q1.py     - quantum solver (QUBO -> Ising Hamiltonian -> QAOA, p=2
            layers, optimized with COBYLA, sampled 1024 shots)

How to run
----------
1. Clone the repo:
   git clone https://github.com/Abdul-Rehman-9127/Pakquant.git
   cd Pakquant

2. Create and activate a virtual environment:
   python3 -m venv .venv
   source .venv/bin/activate

3. Install dependencies:
   pip install qiskit scipy numpy

4. Run the classical baseline:
   python brute.py
   Expected output: best k = 7  score = 220  bits = 00111
   (bits = 00111 means ICU + dialysis + water are kept powered)

5. Run the quantum solver:
   python q1.py
   This prints the Hamiltonian, runs a sanity check against the cost
   function, runs the QAOA optimization loop (takes roughly 1-2
   minutes on a laptop), then prints the top 10 measurement results.
   Expected top result: 11100 (around 29% of 1024 shots), which
   matches the brute-force answer above (Qiskit prints bitstrings in
   reversed order, so 11100 here corresponds to bits 00111 above).

Results summary
----------------
Both the classical baseline and the quantum solver agree on the same
optimal answer: keep ICU, dialysis, and water pumping station powered,
for a total lives-at-risk score of 220. QAOA reaches this answer on
about 29% of individual measurement shots at a shallow circuit depth
(p=2), run on a classical simulator (not real quantum hardware).

Credits
-------
Built using Qiskit (IBM open-source quantum SDK) and SciPy (COBYLA
optimizer). AI assistance (Claude) was used for explanation and
debugging support while writing and understanding this code. All code
was written and understood by the team, and all results are our own.
