#!/usr/bin/env python3
"""
测试 ω=1.0 的收敛行为
"""

import numpy as np
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)
omega = 1.0
wave = WaveParams(omega=omega, l=2)
solver = RadialFlowSolver(bh, wave)

print(f"测试 ω = {omega}")
print("=" * 60)

# 使用更宽松的阈值
result = solver.solve(
    r_max=200.0,
    convergence_threshold=0.2,  # 20% 阈值（宽松）
    n_periods=2,
    max_r_factor=10.0
)

print(f"最终 r/r0 = {result['r'][-1]/bh.r0:.1f}")
print(f"积分点数 = {len(result['r'])}")
print(f"收敛状态 = {result['success']}")
print(f"|A/B| = {np.abs(result['reflection_coeff']):.6f}")
print()

# 分析最后一段
n_check = min(5000, len(result['sigma']))
sigma_last = result['sigma'][-n_check:]
sigma_mean = np.mean(sigma_last)
deviations = np.abs(sigma_last - sigma_mean)
max_deviation = np.max(deviations)
relative_osc = max_deviation / np.abs(sigma_mean) if np.abs(sigma_mean) > 0 else np.inf

print(f"最后 {n_check} 点:")
print(f"  平均值 = {sigma_mean}")
print(f"  |平均值| = {np.abs(sigma_mean):.6f}")
print(f"  最大偏差 = {max_deviation:.6f}")
print(f"  相对震荡 = {relative_osc:.6f} ({relative_osc*100:.2f}%)")
