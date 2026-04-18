import csv
import numpy as np

with open('comparison_results.csv', 'r') as f:
    reader = csv.DictReader(f)
    data = list(reader)

# Separate cases with/without reference
with_ref = [r for r in data if r['ref_ratio']]
without_ref = [r for r in data if not r['ref_ratio']]

print(f"Total cases: {len(data)}")
print(f"Cases with reference: {len(with_ref)}")
print(f"Cases without reference (MST failed): {len(without_ref)}\n")

# Analyze agreement
diffs = [float(r['diff']) for r in with_ref]
refs = [float(r['ref_ratio']) for r in with_ref]
rel_errors = [abs(d/r) if r > 1e-10 else 0 for d, r in zip(diffs, refs)]

print("Agreement statistics (cases with reference):")
print(f"  Max absolute diff: {max(diffs):.3e}")
print(f"  Mean absolute diff: {np.mean(diffs):.3e}")
print(f"  Max relative error: {max(rel_errors):.2%}")
print(f"  Mean relative error: {np.mean(rel_errors):.2%}\n")

# Find worst cases
worst_idx = np.argmax(diffs)
print(f"Worst case (absolute diff):")
r = with_ref[worst_idx]
print(f"  a={r['a']}, ω={r['omega']}, l={r['l']}, m={r['m']}")
print(f"  ref={float(r['ref_ratio']):.6e}, calc={float(r['calc_ratio']):.6e}\n")

# Show MST failed cases
print("Cases where MST method failed (but we computed):")
for r in without_ref:
    print(f"  a={r['a']}, ω={r['omega']}, l={r['l']}, m={r['m']} | calc={float(r['calc_ratio']):.6e}")
