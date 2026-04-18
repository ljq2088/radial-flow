#!/usr/bin/env python3
"""
测试简化版本的相位修正：r vs r*
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
print("简化版本的相位修正对比：r vs r*")
print("=" * 80)
print(f"参数: ω = {omega}, l = {l}, r₀ = {bh.r0}")
print()

# 计算简化版本
wave_simple = WaveParams(omega=omega, l=l, use_subleading_term=False)
solver_simple = RadialFlowSolver(bh, wave_simple)

# 使用 C++ 求解器直接获取结果
sys.path.insert(0, os.path.join(current_dir, 'python'))
import _radial_flow_cpp as cpp
cpp_solver = cpp.RadialFlowSolver(bh._cpp_obj, wave_simple._cpp_obj)
cpp_result = cpp_solver.solve(200.0, 1e-6, 10000, 0.01, 2, 2.0)

# 转换为 numpy 数组
r = np.array(cpp_result.r)
sigma_raw = np.array(cpp_result.sigma)
r0 = bh.r0

print("原始 sigma（未修正）:")
print(f"  最后一个点: σ = {sigma_raw[-1]:.6f}")
print(f"  -σ = {-sigma_raw[-1]:.6f}")
print()

# 方法1: 使用龟坐标 r* 修正（当前代码）
print("方法1: 使用龟坐标 r* 修正")
print("-" * 80)
r_star = r + r0 * np.log(r / r0 - 1.0)
phase_correction_rstar = np.exp(-2.0j * omega * r_star)
sigma_corrected_rstar = sigma_raw * phase_correction_rstar
A_B_rstar = -sigma_corrected_rstar[-1]

print(f"  r* = r + r₀ ln(r/r₀ - 1)")
print(f"  最后一点: r = {r[-1]:.2f}, r* = {r_star[-1]:.2f}")
print(f"  相位修正: exp(-2iω r*) = {phase_correction_rstar[-1]:.6f}")
print(f"  修正后: σ_corrected = {sigma_corrected_rstar[-1]:.6f}")
print(f"  反射系数: A/B = {A_B_rstar:.6f}")
print(f"  |A/B| = {abs(A_B_rstar):.6f}")
print()

# 方法2: 使用普通坐标 r 修正（简化版本应该用这个？）
print("方法2: 使用普通坐标 r 修正")
print("-" * 80)
phase_correction_r = np.exp(-2.0j * omega * r)
sigma_corrected_r = sigma_raw * phase_correction_r
A_B_r = -sigma_corrected_r[-1]

print(f"  直接使用 r")
print(f"  最后一点: r = {r[-1]:.2f}")
print(f"  相位修正: exp(-2iω r) = {phase_correction_r[-1]:.6f}")
print(f"  修正后: σ_corrected = {sigma_corrected_r[-1]:.6f}")
print(f"  反射系数: A/B = {A_B_r:.6f}")
print(f"  |A/B| = {abs(A_B_r):.6f}")
print()

# 方法3: 不做相位修正
print("方法3: 不做相位修正")
print("-" * 80)
A_B_no_correction = -sigma_raw[-1]
print(f"  直接使用原始 σ")
print(f"  反射系数: A/B = {A_B_no_correction:.6f}")
print(f"  |A/B| = {abs(A_B_no_correction):.6f}")
print()

# 对比
print("=" * 80)
print("对比分析")
print("=" * 80)
print(f"{'方法':<30} {'|A/B|':<15} {'相对差异':<15}")
print("-" * 80)

ref = abs(A_B_rstar)
print(f"{'方法1: r* 修正（当前）':<30} {abs(A_B_rstar):<15.6f} {'0.00%':<15}")
print(f"{'方法2: r 修正':<30} {abs(A_B_r):<15.6f} {(abs(A_B_r)-ref)/ref*100:>14.2f}%")
print(f"{'方法3: 无修正':<30} {abs(A_B_no_correction):<15.6f} {(abs(A_B_no_correction)-ref)/ref*100:>14.2f}%")

print()
print("=" * 80)
print("物理解释")
print("=" * 80)
print()
print("简化版本 (V = l(l+1)/r²):")
print("  渐近行为: R ~ A e^{iωr} + B e^{-iωr}")
print("  应该使用普通坐标 r，而不是龟坐标 r*")
print()
print("完整版本 (V = l(l+1)/r² + f'/r):")
print("  渐近行为: R ~ (A₀ + A₁/r) r^{-1} e^{iωr*}")
print("  应该使用龟坐标 r*")
print()

# 计算 r 和 r* 的差异
print("r 和 r* 的差异:")
print(f"  r = {r[-1]:.2f}")
print(f"  r* = {r_star[-1]:.2f}")
print(f"  r* - r = {r_star[-1] - r[-1]:.2f}")
print(f"  相位差: 2ω(r* - r) = {2*omega*(r_star[-1] - r[-1]):.2f} rad")
print(f"  相位差: {2*omega*(r_star[-1] - r[-1])/np.pi:.2f} π")
