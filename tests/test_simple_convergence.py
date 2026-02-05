#!/usr/bin/env python3
"""
简单测试：手动检查收敛判据
"""

import numpy as np
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)
omega = 0.5  # 测试 ω=0.5
wave = WaveParams(omega=omega, l=2)
solver = RadialFlowSolver(bh, wave)

print(f"测试 ω = {omega}")
print("=" * 60)

# 使用更宽松的阈值
result = solver.solve(
    r_max=200.0,
    convergence_threshold=0.2,  # 20% 阈值（非常宽松）
    n_periods=2,
    max_r_factor=10.0
)

print(f"最终 r/r0 = {result['r'][-1]/bh.r0:.1f}")
print(f"积分点数 = {len(result['r'])}")
print(f"收敛状态 = {result['success']}")
print(f"|A/B| = {np.abs(result['reflection_coeff']):.6f}")
