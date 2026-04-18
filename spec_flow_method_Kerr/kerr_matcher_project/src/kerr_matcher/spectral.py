from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .chebyshev import ChebyshevGrid
from .kerr import Dbar_of_z, Pbar_of_z, r_star, rho1_from_sigma
from .params import SolverParams

ArrayC = NDArray[np.complex128]


@dataclass(slots=True)
class SpectralState:
    grid: ChebyshevGrid
    up: ArrayC
    um: ArrayC
    bdry_plus: complex
    ratio_up_over_um: complex
    rho_match: complex



def _diag(v: ArrayC) -> ArrayC:
    return np.diag(np.asarray(v, dtype=complex))



def A_plus(z: ArrayC | float, params: SolverParams) -> ArrayC | complex:
    zz = np.asarray(z, dtype=complex)
    return zz**-3 * np.exp(1.0j * params.omega * r_star(1.0 / zz, params))



def A_minus(z: ArrayC | float, params: SolverParams) -> ArrayC | complex:
    zz = np.asarray(z, dtype=complex)
    return zz * np.exp(-1.0j * params.omega * r_star(1.0 / zz, params))



def A_plus_logder(z: ArrayC | float, params: SolverParams) -> ArrayC | complex:
    zz = np.asarray(z, dtype=complex)
    return -3.0 / zz - 1.0j * params.omega * Pbar_of_z(zz, params) / (Dbar_of_z(zz, params) * zz**2)



def A_minus_logder(z: ArrayC | float, params: SolverParams) -> ArrayC | complex:
    zz = np.asarray(z, dtype=complex)
    return 1.0 / zz + 1.0j * params.omega * Pbar_of_z(zz, params) / (Dbar_of_z(zz, params) * zz**2)



def coeff_B_plus(z: ArrayC, params: SolverParams) -> ArrayC:
    zz = np.asarray(z, dtype=complex)
    return 6.0 * params.M * zz**2 - 2.0j * params.a**2 * params.omega * zz**2 - 4.0 * params.a**2 * zz**3 - 2.0j * params.omega - 2.0 * zz



def coeff_B_minus(z: ArrayC, params: SolverParams) -> ArrayC:
    zz = np.asarray(z, dtype=complex)
    return -10.0 * params.M * zz**2 + 2.0j * params.a**2 * params.omega * zz**2 + 4.0 * params.a**2 * zz**3 + 2.0j * params.omega + 6.0 * zz



def coeff_C_plus(z: ArrayC, params: SolverParams) -> ArrayC:
    zz = np.asarray(z, dtype=complex)
    M = params.M
    a = params.a
    w = params.omega
    m = params.m_mode
    lam = params.lambda_value
    D = Dbar_of_z(zz, params)
    num = (
        12.0 * M**2 * zz**2
        - 12.0j * M * a**2 * w * zz**2
        - 18.0 * M * a**2 * zz**3
        + 4.0j * M * a * m * zz**2
        + 2.0 * M * lam * zz
        - 6.0 * M * zz
        + 6.0j * a**4 * w * zz**3
        + 6.0 * a**4 * zz**4
        - 2.0 * a**3 * m * w * zz**2
        - a**2 * lam * zz**2
        + a**2 * m**2 * zz**2
        + 6.0j * a**2 * w * zz
        + 6.0 * a**2 * zz**2
        - 2.0 * a * m * w
        - 4.0j * a * m * zz
        - lam
    )
    return num / D



def coeff_C_minus(z: ArrayC, params: SolverParams) -> ArrayC:
    zz = np.asarray(z, dtype=complex)
    M = params.M
    a = params.a
    w = params.omega
    m = params.m_mode
    lam = params.lambda_value
    D = Dbar_of_z(zz, params)
    num = (
        12.0 * M**2 * zz**2
        - 12.0j * M * a**2 * w * zz**2
        - 10.0 * M * a**2 * zz**3
        + 4.0j * M * a * m * zz**2
        + 2.0 * M * lam * zz
        + 8.0j * M * w
        - 14.0 * M * zz
        + 2.0j * a**4 * w * zz**3
        + 2.0 * a**4 * zz**4
        - 2.0 * a**3 * m * w * zz**2
        - a**2 * lam * zz**2
        + a**2 * m**2 * zz**2
        + 2.0j * a**2 * w * zz
        + 6.0 * a**2 * zz**2
        - 2.0 * a * m * w
        - 4.0j * a * m * zz
        - lam
        + 4.0
    )
    return num / D



def build_plus_operator(grid: ChebyshevGrid, params: SolverParams) -> ArrayC:
    z = np.asarray(grid.z, dtype=complex)
    Dbar = Dbar_of_z(z, params)
    B = coeff_B_plus(z, params)
    C = coeff_C_plus(z, params)
    return _diag(z**2 * Dbar) @ grid.D2 + _diag(B) @ grid.D + _diag(C)



def build_minus_operator(grid: ChebyshevGrid, params: SolverParams) -> ArrayC:
    z = np.asarray(grid.z, dtype=complex)
    Dbar = Dbar_of_z(z, params)
    B = coeff_B_minus(z, params)
    C = coeff_C_minus(z, params)
    return _diag(z**2 * Dbar) @ grid.D2 + _diag(B) @ grid.D + _diag(C)



def solve_up(grid: ChebyshevGrid, params: SolverParams) -> ArrayC:
    A = build_plus_operator(grid, params)
    rhs = np.zeros(grid.n, dtype=complex)
    rhs[0] = 1.0 + 0.0j
    A[0, :] = 0.0
    A[0, 0] = 1.0
    return np.linalg.solve(A, rhs)



def boundary_plus_value(up: ArrayC, rho1: complex, grid: ChebyshevGrid, params: SolverParams) -> complex:
    z1 = grid.zc
    d1u = grid.D[0, :] @ up
    u1 = up[0]
    val = A_plus(z1, params) * (d1u + (A_plus_logder(z1, params) - rho1) * u1)
    return complex(val)



def boundary_minus_row(rho1: complex, grid: ChebyshevGrid, params: SolverParams) -> ArrayC:
    z1 = grid.zc
    row = A_minus(z1, params) * grid.D[0, :].astype(complex).copy()
    row[0] += A_minus(z1, params) * (A_minus_logder(z1, params) - rho1)
    return row



def solve_um(grid: ChebyshevGrid, rho1: complex, bdry_plus: complex, params: SolverParams) -> ArrayC:
    A = build_minus_operator(grid, params)
    rhs = np.zeros(grid.n, dtype=complex)
    rhs[0] = bdry_plus
    A[0, :] = boundary_minus_row(rho1, grid, params)
    return np.linalg.solve(A, rhs)



def run_spectral_matching(grid: ChebyshevGrid, sigma1: complex, params: SolverParams) -> SpectralState:
    rho1 = rho1_from_sigma(sigma1, params)
    up = solve_up(grid, params)
    bdry = boundary_plus_value(up, rho1, grid, params)
    um = solve_um(grid, rho1, bdry, params)
    ratio = up[-1] / um[-1]
    return SpectralState(grid=grid, up=up, um=um, bdry_plus=bdry, ratio_up_over_um=ratio, rho_match=rho1)
