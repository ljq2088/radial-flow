#!/usr/bin/env python3
"""
测试高频收敛：ω = 3.0, 5.0
尝试不同的 n_periods 和 threshold
"""

import numpy as np
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)

print("=" * 70)
print("高频收敛测试")
print("=" * 70)
print()

# 测试 ω = 3.0
omega = 3.0
wave = WaveParams(omega=omega, l=2)
solver = RadialFlowSolver(bh, wave)

print(f"测试 ω = {omega}")
print("-" * 70)

# 测试不同的 n_periods
for n_periods in [5, 10, 20, 50]:
    print(f"\nn_periods = {n_periods}:")
    result = solver.solve(
        r_max=200.0,
        convergence_threshold=0.1,  # 10% 阈值
        n_periods=n_periods,
        max_r_factor=10.0
    )

    r_final = result['r'][-1] / bh.r0
    n_points = len(result['r'])
    converged = result['success']
    abs_reflection = np.abs(result['reflection_coeff'])

    print(f"  r_final/r0 = {r_final:.1f}")
    print(f"  积分点数 = {n_points}")
    print(f"  收敛状态 = {converged}")
    print(f"  |A/B| = {abs_reflection:.6f}")

print()
print("=" * 70)

# 测试 ω = 5.0
omega = 5.0
wave = WaveParams(omega=omega, l=2)
solver = RadialFlowSolver(bh, wave)

print(f"测试 ω = {omega}")
print("-" * 70)

# 测试不同的 n_periods
for n_periods in [5, 10, 20, 50]:
    print(f"\nn_periods = {n_periods}:")
    result = solver.solve(
        r_max=200.0,
        convergence_threshold=0.1,  # 10% 阈值
        n_periods=n_periods,
        max_r_factor=10.0
    )

    r_final = result['r'][-1] / bh.r0
    n_points = len(result['r'])
    converged = result['success']
    abs_reflection = np.abs(result['reflection_coeff'])

    print(f"  r_final/r0 = {r_final:.1f}")
    print(f"  积分点数 = {n_points}")
    print(f"  收敛状态 = {converged}")
    print(f"  |A/B| = {abs_reflection:.6f}")

print()
print("=" * 70)
print("测试完成")
print("=" * 70)
