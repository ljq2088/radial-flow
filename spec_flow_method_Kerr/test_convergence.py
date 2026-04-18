import sys
sys.path.insert(0, 'kerr_matcher_project/src')

from kerr_matcher.params import SolverParams
from kerr_matcher.solver import solve_case

# Test one worst case
case = {'a': 0.7, 'omega': 0.5, 'l': 2, 'm': 2, 'lambda': 1.6966, 'ref': 0.021204}

print(f"Case: a={case['a']}, ω={case['omega']}, l={case['l']}, m={case['m']}")
print(f"λ={case['lambda']:.4f}, ref={case['ref']:.6e}\n")

configs = [
    {'n_cheb': 32, 'flow_eps': 1e-6, 'r_match': 8.0},
    {'n_cheb': 48, 'flow_eps': 1e-6, 'r_match': 8.0},
    {'n_cheb': 64, 'flow_eps': 1e-6, 'r_match': 8.0},
    {'n_cheb': 32, 'flow_eps': 1e-8, 'r_match': 8.0},
    {'n_cheb': 32, 'flow_eps': 1e-10, 'r_match': 8.0},
    {'n_cheb': 64, 'flow_eps': 1e-8, 'r_match': 8.0},
]

for cfg in configs:
    try:
        params = SolverParams(
            M=1.0, a=case['a'], omega=case['omega'],
            ell=case['l'], m_mode=case['m'],
            lambda_sep=case['lambda'],
            r_match=cfg['r_match'], n_cheb=cfg['n_cheb'], flow_eps=cfg['flow_eps']
        )
        result = solve_case(params)
        calc = abs(result.spectral.ratio_up_over_um)
        err = abs(calc - case['ref']) / case['ref'] * 100
        print(f"n_cheb={cfg['n_cheb']:2}, flow_eps={cfg['flow_eps']:.0e}: {calc:.6e} (err={err:6.2f}%)")
    except Exception as e:
        print(f"n_cheb={cfg['n_cheb']:2}, flow_eps={cfg['flow_eps']:.0e}: FAILED")
