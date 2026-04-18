from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .params import SolverParams


@dataclass(slots=True)
class FlowResult:
    z_start: float
    z_end: float
    Sigma_start: complex
    Sigma_end: complex
    sigma_logderiv: complex
    nfev: int
    success: bool
    message: str


def horizon_sigma_prime(params: SolverParams) -> complex:
    r"""\Sigma'(z_h) 的视界极限。"""
    r0 = params.r0
    w = params.omega
    ell = params.ell
    num = 4.0 * w**2 * r0**4 - ell * (ell + 1.0) * r0**2
    den = -1j * w - 2.0 * w**2 * r0
    return num / den


def sigma_flow_rhs(z: float, Sigma: complex, params: SolverParams) -> complex:
    r0 = params.r0
    w = params.omega
    ell = params.ell
    f = 1.0 - r0 * z
    term1 = (-w**2 * z**-2 / f + ell * (ell + 1.0)) / (1j * w * z**2)
    term2 = -(1j * w / f) * Sigma**2
    return term1 + term2


def sigma_from_Sigma(z: float, Sigma: complex, params: SolverParams) -> complex:
    return 1j * params.omega * Sigma / (1.0 - params.r0 * z)


def integrate_horizon_flow(params: SolverParams) -> FlowResult:
    z_h = params.z_horizon
    if not (0.0 < params.z_match < z_h):
        raise ValueError(f"需要满足 0 < z_match < z_h={z_h}, 当前 z_match={params.z_match}")

    z_start = z_h - params.flow_eps
    dSigma_h = horizon_sigma_prime(params)
    Sigma_start = params.r0**2 + dSigma_h * (z_start - z_h)

    def rhs(z: float, y: np.ndarray) -> np.ndarray:
        Sigma = y[0] + 1j * y[1]
        dSigma = sigma_flow_rhs(z, Sigma, params)
        return np.array([dSigma.real, dSigma.imag], dtype=float)

    sol = solve_ivp(
        rhs,
        t_span=(z_start, params.z_match),
        y0=np.array([Sigma_start.real, Sigma_start.imag], dtype=float),
        method="DOP853",
        rtol=params.flow_rtol,
        atol=params.flow_atol,
    )
    Sigma_end = sol.y[0, -1] + 1j * sol.y[1, -1]
    sigma_log = sigma_from_Sigma(params.z_match, Sigma_end, params)

    return FlowResult(
        z_start=z_start,
        z_end=params.z_match,
        Sigma_start=Sigma_start,
        Sigma_end=Sigma_end,
        sigma_logderiv=sigma_log,
        nfev=sol.nfev,
        success=sol.success,
        message=sol.message,
    )
