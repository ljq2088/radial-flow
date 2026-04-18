from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .params import SolverParams

ArrayC = NDArray[np.complex128]
ArrayF = NDArray[np.float64]


def Delta(r: float | ArrayF | ArrayC, params: SolverParams) -> complex | ArrayF | ArrayC:
    return np.asarray(r) ** 2 - 2.0 * params.M * np.asarray(r) + params.a**2


def Delta_prime(r: float | ArrayF | ArrayC, params: SolverParams) -> complex | ArrayF | ArrayC:
    return 2.0 * np.asarray(r) - 2.0 * params.M


def P_of_r(r: float | ArrayF | ArrayC, params: SolverParams) -> complex | ArrayF | ArrayC:
    return np.asarray(r) ** 2 + params.a**2


def K_of_r(r: float | ArrayF | ArrayC, params: SolverParams) -> complex | ArrayF | ArrayC:
    rr = np.asarray(r)
    return (rr**2 + params.a**2) * params.omega - params.a * params.m_mode


def potential_V(r: float | ArrayF | ArrayC, params: SolverParams) -> complex | ArrayF | ArrayC:
    rr = np.asarray(r)
    D = Delta(rr, params)
    K = K_of_r(rr, params)
    lam = params.lambda_value
    return (K**2 + 4.0j * (rr - params.M) * K) / D - 8.0j * params.omega * rr - lam


def r_star(r: float | ArrayF | ArrayC, params: SolverParams) -> complex | ArrayF | ArrayC:
    rr = np.asarray(r, dtype=complex)
    rp = params.r_plus
    rm = params.r_minus
    coeff_p = 2.0 * params.M * rp / (rp - rm)
    coeff_m = 2.0 * params.M * rm / (rp - rm)
    return rr + coeff_p * np.log(rr - rp) - coeff_m * np.log(rr - rm)


def Dbar_of_z(z: float | ArrayF | ArrayC, params: SolverParams) -> complex | ArrayF | ArrayC:
    zz = np.asarray(z)
    return 1.0 - 2.0 * params.M * zz + params.a**2 * zz**2


def Pbar_of_z(z: float | ArrayF | ArrayC, params: SolverParams) -> complex | ArrayF | ArrayC:
    zz = np.asarray(z)
    return 1.0 + params.a**2 * zz**2


def Kbar_of_z(z: float | ArrayF | ArrayC, params: SolverParams) -> complex | ArrayF | ArrayC:
    zz = np.asarray(z)
    return params.omega * Pbar_of_z(zz, params) - params.a * params.m_mode * zz**2


def Vz_of_z(z: float | ArrayF | ArrayC, params: SolverParams) -> complex | ArrayF | ArrayC:
    zz = np.asarray(z, dtype=complex)
    D = Dbar_of_z(zz, params)
    Kb = Kbar_of_z(zz, params)
    lam = params.lambda_value
    return Kb**2 / (D * zz**2) + 4.0j * (1.0 - params.M * zz) * Kb / (D * zz) - 8.0j * params.omega / zz - lam


def rho1_from_sigma(sigma1: complex, params: SolverParams) -> complex:
    z1 = params.z_match
    Pbar1 = Pbar_of_z(z1, params)
    Dbar1 = Dbar_of_z(z1, params)
    return -1.0j * params.k_horizon * Pbar1 / (Dbar1 * z1**2) * sigma1
