from __future__ import annotations

from schw_matcher.params import SolverParams
from schw_matcher.solver import format_summary, solve_case


def main() -> None:
    params = SolverParams(r0=2.0, omega=1.2, ell=1, z_match=0.25, n_cheb=25, flow_eps=1.0e-6)
    result = solve_case(params)
    print(format_summary(result))


if __name__ == "__main__":
    main()
