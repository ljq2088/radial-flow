#!/usr/bin/env python3
"""
测试高频时增加积分点数（减小步长）
"""

import numpy as np
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)
omega = 3.0
wave = WaveParams(omega=omega, l=2)
solver = RadialFlowSolver(bh, wave)

print(f"测试 ω = {omega}")
print("=" * 60)

# 测试不同的 n_points
for n_points in [10000, 50000, 100000, 200000]:
    print(f"\nn_points = {n_points}:")

    # 计算步长
    r_start = 2.0 + 1e-6 * 2.0
    r_end = 200.0 * 2.0
    dr = (r_end - r_start) / (n_points - 1)
    print(f"  步长 dr = {dr:.6e}")

    result = solver.solve(
        r_max=200.0,
        convergence_threshold=1e10,  # 强制积分到 r_max
        n_periods=2,
        max_r_factor=1.0,
        n_points=n_points
    )

    r = result['r']
    sigma = result['sigma']

    # 检查是否有 nan
    has_nan = np.any(np.isnan(sigma))
    if has_nan:
        # 找到第一个 nan 的位置
        nan_idx = np.where(np.isnan(sigma))[0][0]
        print(f"  第一个 nan 出现在: idx={nan_idx}, r/r0={r[nan_idx]/bh.r0:.2f}")
        print(f"  σ[{nan_idx-1}] = {sigma[nan_idx-1]}")
        print(f"  |σ[{nan_idx-1}]| = {np.abs(sigma[nan_idx-1]):.6e}")
    else:
        print(f"  无 nan！")
        print(f"  σ[-1] = {sigma[-1]}")
        print(f"  |σ[-1]| = {np.abs(sigma[-1]):.6e}")
        print(f"  反射系数 |A/B| = {np.abs(result['reflection_coeff']):.6e}")

print()
print("=" * 60)
