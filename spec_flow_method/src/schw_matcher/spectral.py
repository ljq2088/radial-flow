from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .chebyshev import ChebyshevGrid
from .params import SolverParams

ArrayC = NDArray[np.complex128]


@dataclass(slots=True)
class SpectralState:
    grid: ChebyshevGrid
    fz: ArrayC
    up: ArrayC
    um: ArrayC
    bdry_plus: complex
    ratio_up_over_um: complex


def phase_S(z: ArrayC | float, params: SolverParams) -> ArrayC | complex:
    r0 = params.r0
    w = params.omega
    return w * (1.0 / z + r0 * np.log(1.0 / (r0 * z)))


def _diag(v: ArrayC) -> ArrayC:
    return np.diag(np.asarray(v, dtype=complex))


def build_plus_operator(grid: ChebyshevGrid, params: SolverParams) -> ArrayC:
    z = np.asarray(grid.z, dtype=complex)
    D = grid.D
    r0 = params.r0
    w = params.omega
    ell = params.ell
    fz = 1.0 - r0 * z

    A = D @ _diag(z**2 * fz) @ D
    A += (-2.0j * w) * _diag(1.0 - r0**2 * z**2) @ D
    V = 2.0j * w * r0**2 * z - r0 * z + (w**2 * r0**2 / fz) * (2.0 - r0**2 * z**2) - ell * (ell + 1.0)
    A += _diag(V)
    return A


def build_minus_operator(grid: ChebyshevGrid, params: SolverParams) -> ArrayC:
    z = np.asarray(grid.z, dtype=complex)
    D = grid.D
    r0 = params.r0
    w = params.omega
    ell = params.ell
    fz = 1.0 - r0 * z

    A = D @ _diag(z**2 * fz) @ D
    A += (+2.0j * w) * _diag(1.0 - r0**2 * z**2) @ D
    V = -2.0j * w * r0**2 * z - r0 * z + (w**2 * r0**2 / fz) * (2.0 - r0**2 * z**2) - ell * (ell + 1.0)
    A += _diag(V)
    return A


def solve_up(grid: ChebyshevGrid, params: SolverParams) -> ArrayC:
    A = build_plus_operator(grid, params)
    rhs = np.zeros(grid.n, dtype=complex)
    rhs[0] = 1.0 + 0.0j
    A[0, :] = 0.0
    A[0, 0] = 1.0
    return np.linalg.solve(A, rhs)


def boundary_plus_value(up: ArrayC, sigma_log: complex, grid: ChebyshevGrid, params: SolverParams) -> complex:
    zc = grid.zc
    r0 = params.r0
    w = params.omega
    d1u = grid.D[0, :] @ up
    u1 = up[0]
    val = (zc**2 * d1u + (zc - 1j * w - 1j * w * r0 * zc - zc**2 * sigma_log) * u1) * np.exp(1j * phase_S(zc, params))
    return complex(val)


def boundary_minus_row(sigma_log: complex, grid: ChebyshevGrid, params: SolverParams) -> ArrayC:
    zc = grid.zc
    r0 = params.r0
    w = params.omega
    row = zc**2 * grid.D[0, :].astype(complex).copy()
    row[0] += zc + 1j * w + 1j * w * r0 * zc - zc**2 * sigma_log
    row *= np.exp(-1j * phase_S(zc, params))
    return row


def solve_um(grid: ChebyshevGrid, sigma_log: complex, bdry_plus: complex, params: SolverParams) -> ArrayC:
    A = build_minus_operator(grid, params)
    rhs = np.zeros(grid.n, dtype=complex)
    rhs[0] = bdry_plus
    A[0, :] = boundary_minus_row(sigma_log, grid, params)
    return np.linalg.solve(A, rhs)


def run_spectral_matching(grid: ChebyshevGrid, sigma_log: complex, params: SolverParams) -> SpectralState:
    z = np.asarray(grid.z, dtype=complex)
    fz = 1.0 - params.r0 * z
    up = solve_up(grid, params)
    bdry = boundary_plus_value(up, sigma_log, grid, params)
    um = solve_um(grid, sigma_log, bdry, params)
    ratio = up[-1] / um[-1]
    return SpectralState(grid=grid, fz=fz, up=up, um=um, bdry_plus=bdry, ratio_up_over_um=ratio)
