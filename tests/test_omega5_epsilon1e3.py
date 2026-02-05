#!/usr/bin/env python3
"""
测试 ω = 5.0, epsilon = 1e-3, 增加 n_points
"""

import numpy as np
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)
omega = 5.0
wave = WaveParams(omega=omega, l=2)
solver = RadialFlowSolver(bh, wave)

print(f"测试 ω = {omega}, epsilon = 1e-3")
print("=" * 60)

epsilon = 1e-3

# 测试不同的 n_points
for n_points in [100000, 200000, 400000, 800000]:
    r_start = bh.r0 * (1 + epsilon)
    r_end = 200.0 * bh.r0
    dr = (r_end - r_start) / (n_points - 1)

    print(f"\nn_points = {n_points}, dr = {dr:.6e}")

    result = solver.solve(
        r_max=200.0,
        epsilon=epsilon,
        convergence_threshold=1e10,
        n_periods=2,
        max_r_factor=1.0,
        n_points=n_points
    )

    sigma = result['sigma']
    has_nan = np.any(np.isnan(sigma))

    if has_nan:
        nan_idx = np.where(np.isnan(sigma))[0][0]
        print(f"  第一个 nan 在: idx={nan_idx}, r/r0={result['r'][nan_idx]/bh.r0:.2f}")
    else:
        print(f"  无 nan！|A/B| = {np.abs(result['reflection_coeff']):.6e}")

print()
print("=" * 60)
