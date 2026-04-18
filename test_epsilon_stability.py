#!/usr/bin/env python3
"""
测试不同 epsilon 值对结果的影响
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

# 不同的 epsilon 值
epsilon_values = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]

print("=" * 80)
print("测试不同 epsilon 值对结果的影响")
print("=" * 80)
print(f"参数: ω = {omega}, l = {l}, r₀ = {bh.r0}")
print()

# 测试简化版本
print("【简化版本】V = l(l+1)/r²")
print("-" * 80)
print(f"{'epsilon':<12} {'|A/B|':<15} {'arg(A/B)':<15} {'收敛':<8} {'r_end':<10}")
print("-" * 80)

wave_simple = WaveParams(omega=omega, l=l, use_subleading_term=False)
solver_simple = RadialFlowSolver(bh, wave_simple)

results_simple = []
for eps in epsilon_values:
    result = solver_simple.solve(r_max=200.0, epsilon=eps, n_points=10000)
    A_B = result['reflection_coeff']
    results_simple.append({
        'epsilon': eps,
        'A_B': A_B,
        'abs': abs(A_B),
        'arg': np.angle(A_B),
        'success': result['success'],
        'r_end': result['r'][-1]
    })
    print(f"{eps:<12.1e} {abs(A_B):<15.6f} {np.angle(A_B):<15.6f} "
          f"{'Yes' if result['success'] else 'No':<8} {result['r'][-1]:<10.2f}")

# 计算简化版本的标准差
abs_values_simple = [r['abs'] for r in results_simple]
mean_simple = np.mean(abs_values_simple)
std_simple = np.std(abs_values_simple)
rel_std_simple = std_simple / mean_simple * 100

print()
print(f"统计分析:")
print(f"  平均值: {mean_simple:.6f}")
print(f"  标准差: {std_simple:.6e}")
print(f"  相对标准差: {rel_std_simple:.4f}%")
print(f"  最大值: {max(abs_values_simple):.6f}")
print(f"  最小值: {min(abs_values_simple):.6f}")
print(f"  极差: {max(abs_values_simple) - min(abs_values_simple):.6e}")

print()
print()

# 测试完整版本
print("【完整版本】V = l(l+1)/r² + f'/r")
print("-" * 80)
print(f"{'epsilon':<12} {'|A/B|':<15} {'arg(A/B)':<15} {'收敛':<8} {'r_end':<10}")
print("-" * 80)

wave_full = WaveParams(omega=omega, l=l, use_subleading_term=True)
solver_full = RadialFlowSolver(bh, wave_full)

results_full = []
for eps in epsilon_values:
    result = solver_full.solve(r_max=200.0, epsilon=eps, n_points=10000)
    A_B = result['reflection_coeff']
    results_full.append({
        'epsilon': eps,
        'A_B': A_B,
        'abs': abs(A_B),
        'arg': np.angle(A_B),
        'success': result['success'],
        'r_end': result['r'][-1]
    })
    print(f"{eps:<12.1e} {abs(A_B):<15.6f} {np.angle(A_B):<15.6f} "
          f"{'Yes' if result['success'] else 'No':<8} {result['r'][-1]:<10.2f}")

# 计算完整版本的标准差
abs_values_full = [r['abs'] for r in results_full]
mean_full = np.mean(abs_values_full)
std_full = np.std(abs_values_full)
rel_std_full = std_full / mean_full * 100

print()
print(f"统计分析:")
print(f"  平均值: {mean_full:.6f}")
print(f"  标准差: {std_full:.6e}")
print(f"  相对标准差: {rel_std_full:.4f}%")
print(f"  最大值: {max(abs_values_full):.6f}")
print(f"  最小值: {min(abs_values_full):.6f}")
print(f"  极差: {max(abs_values_full) - min(abs_values_full):.6e}")

print()
print("=" * 80)
print("结论")
print("=" * 80)

if rel_std_simple < 0.1 and rel_std_full < 0.1:
    print("✅ 结果非常稳定！不同 epsilon 值给出一致的结果（相对标准差 < 0.1%）")
elif rel_std_simple < 1.0 and rel_std_full < 1.0:
    print("✅ 结果稳定。不同 epsilon 值给出基本一致的结果（相对标准差 < 1%）")
else:
    print("⚠️  结果对 epsilon 值敏感，可能需要调整数值方法")

print()
print(f"简化版本相对标准差: {rel_std_simple:.4f}%")
print(f"完整版本相对标准差: {rel_std_full:.4f}%")
