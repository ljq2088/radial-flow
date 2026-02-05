#!/usr/bin/env python3
"""
测试 ω=1.0 使用更多周期
"""

import numpy as np
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)
omega = 1.0
wave = WaveParams(omega=omega, l=2)
solver = RadialFlowSolver(bh, wave)

print(f"测试 ω = {omega}")
print("=" * 60)

# 测试不同的 n_periods
for n_periods in [2, 5, 10]:
    print(f"\nn_periods = {n_periods}:")
    result = solver.solve(
        r_max=200.0,
        convergence_threshold=0.1,  # 10% 阈值
        n_periods=n_periods,
        max_r_factor=10.0
    )

    print(f"  最终 r/r0 = {result['r'][-1]/bh.r0:.1f}")
    print(f"  积分点数 = {len(result['r'])}")
    print(f"  收敛状态 = {result['success']}")
    print(f"  |A/B| = {np.abs(result['reflection_coeff']):.6f}")
