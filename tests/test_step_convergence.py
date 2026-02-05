#!/usr/bin/env python3
"""
测试高频时的步长收敛性
"""

import numpy as np
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)

print("=" * 70)
print("高频步长收敛性测试")
print("=" * 70)
print()

for omega in [3.0, 5.0]:
    wave = WaveParams(omega=omega, l=2)
    solver = RadialFlowSolver(bh, wave)

    print(f"ω = {omega}")
    print("-" * 70)

    results = []

    # 测试不同的 n_points
    for n_points in [50000, 100000, 200000, 400000]:
        r_start = 2.0 + 1e-6 * 2.0
        r_end = 200.0 * 2.0
        dr = (r_end - r_start) / (n_points - 1)

        result = solver.solve(
            r_max=200.0,
            convergence_threshold=1e10,  # 强制积分到 r_max
            n_periods=2,
            max_r_factor=1.0,
            n_points=n_points
        )

        sigma = result['sigma']
        has_nan = np.any(np.isnan(sigma))

        if not has_nan:
            abs_reflection = np.abs(result['reflection_coeff'])
            results.append({
                'n_points': n_points,
                'dr': dr,
                'abs_reflection': abs_reflection
            })
            print(f"  n_points = {n_points:6d}, dr = {dr:.6e}, |A/B| = {abs_reflection:.6f}")
        else:
            print(f"  n_points = {n_points:6d}, dr = {dr:.6e}, 包含 nan")

    # 检查收敛性
    if len(results) >= 2:
        print(f"\n  收敛性检查:")
        for i in range(1, len(results)):
            diff = abs(results[i]['abs_reflection'] - results[i-1]['abs_reflection'])
            rel_diff = diff / results[i]['abs_reflection'] if results[i]['abs_reflection'] > 0 else 0
            print(f"    Δ|A/B| ({results[i-1]['n_points']} → {results[i]['n_points']}) = {diff:.6e} ({rel_diff*100:.2f}%)")

    print()

print("=" * 70)
