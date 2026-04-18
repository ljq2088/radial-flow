import sys
import csv
import numpy as np
sys.path.insert(0, 'kerr_matcher_project/src')

from kerr_matcher.params import SolverParams
from kerr_matcher.solver import solve_case

# Read CSV
with open('amplitude_ratio_table.csv', 'r') as f:
    reader = csv.DictReader(f)
    cases = list(reader)

results = []
for i, row in enumerate(cases):
    a = float(row['a'])
    omega = float(row['omega'])
    ell = int(row['l'])
    m = int(row['m'])
    lam = float(row['lambda'])

    # Handle missing reference data
    ref_ratio = float(row['|B_ref/B_inc|']) if row['|B_ref/B_inc|'] else None

    try:
        params = SolverParams(
            M=1.0, a=a, omega=omega, ell=ell, m_mode=m,
            lambda_sep=lam, r_match=8.0, n_cheb=32, flow_eps=1e-6
        )
        result = solve_case(params)
        calc_ratio = abs(result.spectral.ratio_up_over_um)

        diff = abs(calc_ratio - ref_ratio) if ref_ratio is not None else None
        results.append({
            'a': a, 'omega': omega, 'l': ell, 'm': m, 'lambda': lam,
            'ref_ratio': ref_ratio if ref_ratio else '',
            'calc_ratio': calc_ratio,
            'diff': diff if diff is not None else ''
        })

        if ref_ratio:
            print(f"Row {i+2}: a={a}, ω={omega}, l={ell}, m={m} | ref={ref_ratio:.6e}, calc={calc_ratio:.6e}")
        else:
            print(f"Row {i+2}: a={a}, ω={omega}, l={ell}, m={m} | ref=N/A, calc={calc_ratio:.6e}")
    except Exception as e:
        print(f"Row {i+2} failed: {e}")

# Save results
with open('comparison_results.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['a','omega','l','m','lambda','ref_ratio','calc_ratio','diff'])
    writer.writeheader()
    writer.writerows(results)

print(f"\nCompleted {len(results)} cases. Results saved to comparison_results.csv")
