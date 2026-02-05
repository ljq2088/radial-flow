#!/usr/bin/env python3
"""
测试不同的 epsilon（初始点位置）
"""

import numpy as np
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)

print("=" * 70)
print("测试不同的 epsilon")
print("=" * 70)
print()

for omega in [3.0, 5.0]:
    wave = WaveParams(omega=omega, l=2)
    solver = RadialFlowSolver(bh, wave)

    print(f"ω = {omega}")
    print("-" * 70)

    # 测试不同的 epsilon
    for epsilon in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]:
        r_start = bh.r0 * (1 + epsilon)

        # 计算初值
        l = wave.l
        f0_prime = 1.0 / bh.r0
        numerator = -f0_prime * (1.0 + (l*(l+1) + bh.r0*f0_prime)/(2.0j*omega*bh.r0))
        denominator = bh.r0*f0_prime - 2.0j*omega*bh.r0
        sigma0_prime = numerator / denominator
        sigma_start = 0.0 + sigma0_prime * epsilon * bh.r0

        print(f"\n  epsilon = {epsilon:.0e}:")
        print(f"    r_start/r0 = {r_start/bh.r0:.6f}")
        print(f"    |σ_start| = {np.abs(sigma_start):.6e}")

        result = solver.solve(
            r_max=200.0,
            epsilon=epsilon,
            convergence_threshold=1e10,
            n_periods=2,
            max_r_factor=1.0,
            n_points=100000  # 使用固定的点数
        )

        sigma = result['sigma']
        has_nan = np.any(np.isnan(sigma))

        if has_nan:
            nan_idx = np.where(np.isnan(sigma))[0][0]
            print(f"    第一个 nan 在: idx={nan_idx}, r/r0={result['r'][nan_idx]/bh.r0:.2f}")
        else:
            print(f"    无 nan！|A/B| = {np.abs(result['reflection_coeff']):.6e}")

    print()

print("=" * 70)
