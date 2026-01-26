"""
径向流计算包

C++ 核心 + Python 接口
"""

from .radial_flow import (
    SchwarzschildParams,
    WaveParams,
    RadialFlowSolver,
    plot_sigma_evolution,
    main
)

__version__ = "0.1.0"
__all__ = [
    "SchwarzschildParams",
    "WaveParams",
    "RadialFlowSolver",
    "plot_sigma_evolution",
    "main"
]