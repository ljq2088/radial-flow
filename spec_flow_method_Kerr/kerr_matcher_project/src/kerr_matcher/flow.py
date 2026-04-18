from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .kerr import Delta, Delta_prime, P_of_r, potential_V
from .params import SolverParams


@dataclass(slots=True)
class FlowResult:
    r_start: float
    r_end: float
    sigma_start: complex
    sigma_end: complex
    Rr_over_R_at_match: complex
    nfev: int
    success: bool
    message: str



def sigma_horizon_seed(r: float, params: SolverParams) -> complex:
    """用纯 ingoing 视界主因子 R ~ Δ^2 e^{-ik_H r_*} 给出起始种子。"""
    Dp = Delta_prime(r, params)
    P = P_of_r(r, params)
    return 2.0 * Dp / (1.0j * params.k_horizon * P) - 1.0



def sigma_flow_rhs(r: float, sigma: complex, params: SolverParams) -> complex:
    D = Delta(r, params)
    Dp = Delta_prime(r, params)
    P = P_of_r(r, params)
    V = potential_V(r, params)
    return (2.0 * Dp / D - 2.0 * r / P) * sigma - 1.0j * params.k_horizon * (P / D) * sigma**2 - V / (1.0j * params.k_horizon * P)



def Rr_over_R_from_sigma(r: float, sigma: complex, params: SolverParams) -> complex:
    return 1.0j * params.k_horizon * P_of_r(r, params) / Delta(r, params) * sigma



def integrate_inner_flow(params: SolverParams) -> FlowResult:
    r_start = params.r_start
    sigma_start = sigma_horizon_seed(r_start, params)

    def rhs(r: float, y: np.ndarray) -> np.ndarray:
        sigma = y[0] + 1j * y[1]
        ds = sigma_flow_rhs(r, sigma, params)
        return np.array([ds.real, ds.imag], dtype=float)

    sol = solve_ivp(
        rhs,
        t_span=(r_start, params.r_match),
        y0=np.array([sigma_start.real, sigma_start.imag], dtype=float),
        method="DOP853",
        rtol=params.flow_rtol,
        atol=params.flow_atol,
    )

    sigma_end = sol.y[0, -1] + 1j * sol.y[1, -1]
    logr = Rr_over_R_from_sigma(params.r_match, sigma_end, params)
    return FlowResult(
        r_start=r_start,
        r_end=params.r_match,
        sigma_start=sigma_start,
        sigma_end=sigma_end,
        Rr_over_R_at_match=logr,
        nfev=sol.nfev,
        success=sol.success,
        message=sol.message,
    )
