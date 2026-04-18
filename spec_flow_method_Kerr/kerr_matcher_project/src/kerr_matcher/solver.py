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
        "=== Kerr s=-2 Teukolsky Riccati flow + spectral matching ===",
        f"M = {p.M}",
        f"a = {p.a}",
        f"omega = {p.omega}",
        f"ell = {p.ell}",
        f"m_mode = {p.m_mode}",
        f"lambda = {p.lambda_value}",
        f"r_plus = {p.r_plus}",
        f"Omega_H = {p.omega_horizon}",
        f"k_H = {p.k_horizon}",
        f"r_match = {p.r_match}   (z_match = {p.z_match})",
        f"n_cheb = {p.n_cheb}",
        f"flow start r = {f.r_start}",
        f"flow end   r = {f.r_end}",
        f"sigma(r_match) = {f.sigma_end}",
        f"(R_r/R)(r_match) = {f.Rr_over_R_at_match}",
        f"rho_match = (R_z/R)(z_match) = {result.spectral.rho_match}",
        f"nfev = {f.nfev}",
        f"ratio up/um = {ratio}",
        f"|ratio| = {result.amplitude_ratio_abs}",
        f"arg(ratio) = {result.amplitude_ratio_arg}",
    ]
    return "\n".join(lines)
