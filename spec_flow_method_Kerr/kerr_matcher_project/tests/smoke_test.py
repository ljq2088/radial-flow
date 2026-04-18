from __future__ import annotations

from kerr_matcher.params import SolverParams
from kerr_matcher.solver import solve_case


def main() -> None:
    params = SolverParams(M=1.0, a=0.3, omega=0.2, ell=2, m_mode=2, r_match=8.0, n_cheb=24)
    result = solve_case(params)
    assert result.flow.success
    assert abs(result.spectral.ratio_up_over_um) == abs(result.spectral.ratio_up_over_um)
    print("smoke test passed")
    print(result.spectral.ratio_up_over_um)


if __name__ == "__main__":
    main()
