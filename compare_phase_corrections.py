#!/usr/bin/env python3
"""
对比简化版本和完整版本的相位修正
"""

import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver
import numpy as np

# 测试参数
bh = SchwarzschildParams(r0=2.0)
omega = 0.2
l = 0

print("=" * 80)
print("相位修正对比：简化版本 vs 完整版本")
print("=" * 80)
print(f"参数: ω = {omega}, l = {l}, r₀ = {bh.r0}")
print()

# ============================================================================
# 简化版本
# ============================================================================
print("【简化版本】V = l(l+1)/r²")
print("=" * 80)

wave_simple = WaveParams(omega=omega, l=l, use_subleading_term=False)
solver_simple = RadialFlowSolver(bh, wave_simple)
result_simple = solver_simple._cpp_solver.solve(200.0, 1e-6, 10000, 0.01, 2, 2.0)

r_simple = np.array(result_simple.r)
sigma_simple_raw = np.array(result_simple.sigma)
r0 = bh.r0

# 计算不同的相位修正
r_star_simple = r_simple + r0 * np.log(r_simple / r0 - 1.0)

# r* 修正
phase_rstar_simple = np.exp(-2.0j * omega * r_star_simple[-1])
A_B_rstar_simple = -sigma_simple_raw[-1] * phase_rstar_simple

# r 修正
phase_r_simple = np.exp(-2.0j * omega * r_simple[-1])
A_B_r_simple = -sigma_simple_raw[-1] * phase_r_simple

# 无修正
A_B_no_simple = -sigma_simple_raw[-1]

print(f"r* 修正:  A/B = {A_B_rstar_simple:.6f}, |A/B| = {abs(A_B_rstar_simple):.6f}, arg = {np.angle(A_B_rstar_simple):.6f}")
print(f"r 修正:   A/B = {A_B_r_simple:.6f}, |A/B| = {abs(A_B_r_simple):.6f}, arg = {np.angle(A_B_r_simple):.6f}")
print(f"无修正:   A/B = {A_B_no_simple:.6f}, |A/B| = {abs(A_B_no_simple):.6f}, arg = {np.angle(A_B_no_simple):.6f}")
print()

# ============================================================================
# 完整版本
# ============================================================================
print("【完整版本】V = l(l+1)/r² + f'/r")
print("=" * 80)

wave_full = WaveParams(omega=omega, l=l, use_subleading_term=True)
solver_full = RadialFlowSolver(bh, wave_full)
result_full = solver_full._cpp_solver.solve(200.0, 1e-6, 10000, 0.01, 2, 2.0)

r_full = np.array(result_full.r)
sigma_full_raw = np.array(result_full.sigma)

# 计算不同的相位修正
r_star_full = r_full + r0 * np.log(r_full / r0 - 1.0)

# r* 修正
phase_rstar_full = np.exp(-2.0j * omega * r_star_full[-1])
A_B_rstar_full = -sigma_full_raw[-1] * phase_rstar_full

# r 修正
phase_r_full = np.exp(-2.0j * omega * r_full[-1])
A_B_r_full = -sigma_full_raw[-1] * phase_r_full

# 无修正
A_B_no_full = -sigma_full_raw[-1]

print(f"r* 修正:  A/B = {A_B_rstar_full:.6f}, |A/B| = {abs(A_B_rstar_full):.6f}, arg = {np.angle(A_B_rstar_full):.6f}")
print(f"r 修正:   A/B = {A_B_r_full:.6f}, |A/B| = {abs(A_B_r_full):.6f}, arg = {np.angle(A_B_r_full):.6f}")
print(f"无修正:   A/B = {A_B_no_full:.6f}, |A/B| = {abs(A_B_no_full):.6f}, arg = {np.angle(A_B_no_full):.6f}")
print()

# ============================================================================
# 总结
# ============================================================================
print("=" * 80)
print("总结")
print("=" * 80)
print()

print("关键发现:")
print("  1. 对于 |A/B|（反射系数的模），相位修正方法不影响结果")
print("  2. 对于 arg(A/B)（反射系数的相位），不同修正方法给出不同结果")
print()

print("简化版本:")
print(f"  |A/B| = {abs(A_B_rstar_simple):.6f} (所有方法相同)")
print(f"  arg(A/B): r* 修正 = {np.angle(A_B_rstar_simple):.6f} rad")
print(f"  arg(A/B): r 修正  = {np.angle(A_B_r_simple):.6f} rad")
print(f"  相位差 = {np.angle(A_B_rstar_simple) - np.angle(A_B_r_simple):.6f} rad")
print()

print("完整版本:")
print(f"  |A/B| = {abs(A_B_rstar_full):.6f} (所有方法相同)")
print(f"  arg(A/B): r* 修正 = {np.angle(A_B_rstar_full):.6f} rad")
print(f"  arg(A/B): r 修正  = {np.angle(A_B_r_full):.6f} rad")
print(f"  相位差 = {np.angle(A_B_rstar_full) - np.angle(A_B_r_full):.6f} rad")
print()

print("物理解释:")
print("  - 简化版本理论上应该使用 r 修正（平坦时空渐近）")
print("  - 完整版本应该使用 r* 修正（弯曲时空渐近）")
print("  - 但对于 |A/B|，两种方法给出相同结果")
print("  - 相位差来自 r 和 r* 的差异")
print()

print("结论:")
print("  ✅ 如果只关心 |A/B|，当前实现是正确的")
print("  ⚠️  如果关心 arg(A/B)，简化版本应该使用 r 而不是 r*")
