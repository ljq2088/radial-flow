#!/usr/bin/env python3
"""
测试收敛判据功能
"""

import numpy as np
import matplotlib.pyplot as plt
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

# 测试不同频率
bh = SchwarzschildParams(r0=2.0)
test_omegas = [0.3, 0.5, 1.0, 2.0]

print("=" * 70)
print("测试收敛判据功能")
print("=" * 70)
print()

for omega in test_omegas:
    wave = WaveParams(omega=omega, l=2)
    solver = RadialFlowSolver(bh, wave)

    # 使用收敛判据
    result = solver.solve(
        r_max=200.0,
        convergence_threshold=0.1,  # 1% 阈值
        n_periods=2,  # 检测2个周期
        max_r_factor=5.0  # 最多扩展到 5*r_max
    )

    r = result['r']
    r_final = r[-1] / bh.r0
    n_points = len(r)

    print(f"ω = {omega:.1f}:")
    print(f"  最终 r/r0 = {r_final:.1f}")
    print(f"  积分点数 = {n_points}")
    print(f"  收敛状态 = {result['success']}")
    print(f"  反射系数 A/B = {result['reflection_coeff']:.4f}")
    print(f"  |A/B| = {np.abs(result['reflection_coeff']):.6f}")
    print()

print("=" * 70)
print("测试完成")
print("=" * 70)
