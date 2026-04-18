from __future__ import annotations

from .matching import MatchingResult, run_matching
from .params import SolverParams


def solve_case(params: SolverParams) -> MatchingResult:
    return run_matching(params)


def format_summary(result: MatchingResult) -> str:
    p = result.params
    f = result.flow
    ratio = result.spectral.ratio_up_over_um
    lines = [
        "=== Schwarzschild horizon-flow + spectral matching ===",
        f"r0 = {p.r0}",
        f"omega = {p.omega}",
        f"ell = {p.ell}",
        f"z_match = {p.z_match}   (r_match = {p.r_match})",
        f"n_cheb = {p.n_cheb}",
        f"flow start z = {f.z_start}",
        f"flow end   z = {f.z_end}",
        f"Sigma(z_match) = {f.Sigma_end}",
        f"sigma_log(z_match) = {f.sigma_logderiv}",
        f"nfev = {f.nfev}",
        f"ratio up/um = {ratio}",
        f"|ratio| = {result.amplitude_ratio_abs}",
        f"arg(ratio) = {result.amplitude_ratio_arg}",
    ]
    return "\n".join(lines)
