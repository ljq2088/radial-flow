from __future__ import annotations

import cmath
from dataclasses import dataclass

from .chebyshev import ChebyshevGrid, chebyshev_lobatto_grid
from .flow import FlowResult, integrate_inner_flow
from .params import SolverParams
from .spectral import SpectralState, run_spectral_matching


@dataclass(slots=True)
class MatchingResult:
    params: SolverParams
    flow: FlowResult
    grid: ChebyshevGrid
    spectral: SpectralState

    @property
    def amplitude_ratio_abs(self) -> float:
        return abs(self.spectral.ratio_up_over_um)

    @property
    def amplitude_ratio_arg(self) -> float:
        return cmath.phase(self.spectral.ratio_up_over_um)



def run_matching(params: SolverParams) -> MatchingResult:
    flow = integrate_inner_flow(params)
    if not flow.success:
        raise RuntimeError(f"inner flow 失败: {flow.message}")

    grid = chebyshev_lobatto_grid(params.n_cheb, params.z_match)
    spectral = run_spectral_matching(grid, flow.sigma_end, params)
    return MatchingResult(params=params, flow=flow, grid=grid, spectral=spectral)
