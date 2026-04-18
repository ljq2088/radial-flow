import sys
import csv
sys.path.insert(0, 'kerr_matcher_project/src')

from kerr_matcher.params import SolverParams
from kerr_matcher.solver import solve_case

# Select problematic cases (l=2, m=2, high omega)
test_cases = [
    {'a': 0.0, 'omega': 0.5, 'l': 2, 'm': 2, 'lambda': 4.0, 'ref': 0.024321},
    {'a': 0.3, 'omega': 0.5, 'l': 2, 'm': 2, 'lambda': 3.0058, 'ref': 0.046708},
    {'a': 0.5, 'omega': 0.5, 'l': 2, 'm': 2, 'lambda': 2.3490, 'ref': None},
    {'a': 0.7, 'omega': 0.5, 'l': 2, 'm': 2, 'lambda': 1.6966, 'ref': 0.021204},
    {'a': 0.9, 'omega': 0.5, 'l': 2, 'm': 2, 'lambda': 1.0484, 'ref': 0.036414},
    {'a': 0.9, 'omega': 0.4, 'l': 2, 'm': 2, 'lambda': 1.6316, 'ref': 0.014172},
]

r_match_values = [6.0, 8.0, 10.0, 12.0, 15.0, 20.0]

results = []
for case in test_cases:
    print(f"\nCase: a={case['a']}, ω={case['omega']}, l={case['l']}, m={case['m']}")
    print(f"  λ={case['lambda']:.4f}, ref={case['ref'] if case['ref'] else 'N/A'}")

    for r_match in r_match_values:
        try:
            params = SolverParams(
                M=1.0, a=case['a'], omega=case['omega'],
                ell=case['l'], m_mode=case['m'],
                lambda_sep=case['lambda'],
                r_match=r_match, n_cheb=32, flow_eps=1e-6
            )
            result = solve_case(params)
            calc_ratio = abs(result.spectral.ratio_up_over_um)

            rel_err = None
            if case['ref']:
                rel_err = abs(calc_ratio - case['ref']) / case['ref'] * 100

            results.append({
                'a': case['a'], 'omega': case['omega'],
                'l': case['l'], 'm': case['m'],
                'lambda': case['lambda'], 'ref': case['ref'],
                'r_match': r_match, 'calc_ratio': calc_ratio,
                'rel_err': rel_err
            })

            if rel_err:
                print(f"  r_match={r_match:5.1f}: {calc_ratio:.6e} (err={rel_err:6.2f}%)")
            else:
                print(f"  r_match={r_match:5.1f}: {calc_ratio:.6e}")
        except Exception as e:
            print(f"  r_match={r_match:5.1f}: FAILED - {e}")

# Save results
with open('r_match_scan.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['a','omega','l','m','lambda','ref','r_match','calc_ratio','rel_err'])
    writer.writeheader()
    writer.writerows(results)

print(f"\n\nResults saved to r_match_scan.csv")
