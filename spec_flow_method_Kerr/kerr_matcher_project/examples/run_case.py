from __future__ import annotations

from kerr_matcher.params import SolverParams
from kerr_matcher.solver import format_summary, solve_case


def main() -> None:
    params = SolverParams(
        M=1.0,
        a=0.3,
        omega=0.2,
        ell=2,
        m_mode=2,
        # 如需更精确，请填入实际 spheroidal separation constant
        lambda_sep=3.601,
        r_match=8.0,
        n_cheb=32,
        flow_eps=1.0e-6,
    )
    result = solve_case(params)
    print(format_summary(result))


if __name__ == "__main__":
    main()
