"""PakQuant - classical exact baseline.

Enumerates all 2^5 = 32 on/off combinations of the 5 institutions, keeps the
ones with exactly 3 powered, and reports the best by total lives-at-risk score.

Bit convention (used everywhere in this repo):
    bit i of k  <->  institution i   (i = 0 is ICU ... i = 4 is mall)
    So k = 7 = 0b00111  ->  bits 0,1,2 set  ->  ICU + dialysis + water.
"""

names = ["ICU", "dialysis", "water", "school", "mall"]
scores = [100, 80, 40, 5, 1]
n, k_target = 5, 3

valid = []  # (total_score, k)
for k in range(2 ** n):
    chosen = [i for i in range(n) if (k >> i) & 1]
    if len(chosen) == k_target:
        valid.append((sum(scores[i] for i in chosen), k))

valid.sort(reverse=True)
best_total, best_k = valid[0]

print(f"Valid combinations (exactly {k_target} of {n}): {len(valid)}")
for total, k in valid:
    chosen = [names[i] for i in range(n) if (k >> i) & 1]
    print(f"  k={k:2d}  bits={k:05b}  score={total:3d}  {' + '.join(chosen)}")

print()
print(f"best k = {best_k}  score = {best_total}  bits = {best_k:05b}")
print("powered:", " + ".join(names[i] for i in range(n) if (best_k >> i) & 1))
