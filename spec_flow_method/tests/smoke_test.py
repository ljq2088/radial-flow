from __future__ import annotations

from schw_matcher import SolverParams, solve_case


def main() -> None:
    params = SolverParams(r0=2.0, omega=0.2, ell=0, z_match=0.25, n_cheb=25)
    result = solve_case(params)
    assert result.flow.success
    assert result.grid.n == 25
    assert result.amplitude_ratio_abs >= 0.0
    print("smoke test passed")
    print(result.spectral.ratio_up_over_um)


if __name__ == "__main__":
    main()
