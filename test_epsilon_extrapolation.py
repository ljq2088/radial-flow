#!/usr/bin/env python3
"""
测试更大的 epsilon 值，探索外推到 epsilon=0 的行为
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

# 扩展的 epsilon 值范围（从大到小）
epsilon_values = [1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 1e-6, 1e-7, 1e-8]

print("=" * 80)
print("探索 epsilon → 0 的外推行为")
print("=" * 80)
print(f"参数: ω = {omega}, l = {l}, r₀ = {bh.r0}")
print()

# 测试简化版本
print("【简化版本】V = l(l+1)/r²")
print("-" * 80)
print(f"{'epsilon':<12} {'|A/B|':<15} {'偏差(%)':<12} {'收敛':<8} {'r_end':<10}")
print("-" * 80)

wave_simple = WaveParams(omega=omega, l=l, use_subleading_term=False)
solver_simple = RadialFlowSolver(bh, wave_simple)

results_simple = []
for eps in epsilon_values:
    try:
        result = solver_simple.solve(r_max=200.0, epsilon=eps, n_points=10000)
        A_B = result['reflection_coeff']
        results_simple.append({
            'epsilon': eps,
            'abs': abs(A_B),
            'success': result['success'],
            'r_end': result['r'][-1]
        })
    except Exception as e:
        print(f"{eps:<12.1e} ERROR: {str(e)[:50]}")
        results_simple.append({
            'epsilon': eps,
            'abs': None,
            'success': False,
            'r_end': None
        })

# 计算参考值（使用小 epsilon 的平均值）
ref_simple = np.mean([r['abs'] for r in results_simple[-5:] if r['abs'] is not None])

for r in results_simple:
    if r['abs'] is not None:
        dev = (r['abs'] - ref_simple) / ref_simple * 100
        print(f"{r['epsilon']:<12.1e} {r['abs']:<15.6f} {dev:<12.4f} "
              f"{'Yes' if r['success'] else 'No':<8} {r['r_end']:<10.2f}")
    else:
        print(f"{r['epsilon']:<12.1e} {'FAILED':<15} {'-':<12} {'No':<8} {'-':<10}")

print()
print(f"参考值（小 epsilon 平均）: {ref_simple:.6f}")

print()
print()

# 测试完整版本
print("【完整版本】V = l(l+1)/r² + f'/r")
print("-" * 80)
print(f"{'epsilon':<12} {'|A/B|':<15} {'偏差(%)':<12} {'收敛':<8} {'r_end':<10}")
print("-" * 80)

wave_full = WaveParams(omega=omega, l=l, use_subleading_term=True)
solver_full = RadialFlowSolver(bh, wave_full)

results_full = []
for eps in epsilon_values:
    try:
        result = solver_full.solve(r_max=200.0, epsilon=eps, n_points=10000)
        A_B = result['reflection_coeff']
        results_full.append({
            'epsilon': eps,
            'abs': abs(A_B),
            'success': result['success'],
            'r_end': result['r'][-1]
        })
    except Exception as e:
        print(f"{eps:<12.1e} ERROR: {str(e)[:50]}")
        results_full.append({
            'epsilon': eps,
            'abs': None,
            'success': False,
            'r_end': None
        })

# 计算参考值（使用小 epsilon 的平均值）
ref_full = np.mean([r['abs'] for r in results_full[-5:] if r['abs'] is not None])

for r in results_full:
    if r['abs'] is not None:
        dev = (r['abs'] - ref_full) / ref_full * 100
        print(f"{r['epsilon']:<12.1e} {r['abs']:<15.6f} {dev:<12.4f} "
              f"{'Yes' if r['success'] else 'No':<8} {r['r_end']:<10.2f}")
    else:
        print(f"{r['epsilon']:<12.1e} {'FAILED':<15} {'-':<12} {'No':<8} {'-':<10}")

print()
print(f"参考值（小 epsilon 平均）: {ref_full:.6f}")

print()
print("=" * 80)
print("分析")
print("=" * 80)

# 分析趋势
valid_simple = [(r['epsilon'], r['abs']) for r in results_simple if r['abs'] is not None]
valid_full = [(r['epsilon'], r['abs']) for r in results_full if r['abs'] is not None]

if len(valid_simple) > 1:
    eps_simple, abs_simple = zip(*valid_simple)
    trend_simple = abs_simple[0] - abs_simple[-1]
    print(f"\n简化版本:")
    print(f"  最大 epsilon (ε={eps_simple[0]:.1e}): |A/B| = {abs_simple[0]:.6f}")
    print(f"  最小 epsilon (ε={eps_simple[-1]:.1e}): |A/B| = {abs_simple[-1]:.6f}")
    print(f"  总变化: {trend_simple:.6e} ({trend_simple/abs_simple[-1]*100:.4f}%)")

if len(valid_full) > 1:
    eps_full, abs_full = zip(*valid_full)
    trend_full = abs_full[0] - abs_full[-1]
    print(f"\n完整版本:")
    print(f"  最大 epsilon (ε={eps_full[0]:.1e}): |A/B| = {abs_full[0]:.6f}")
    print(f"  最小 epsilon (ε={eps_full[-1]:.1e}): |A/B| = {abs_full[-1]:.6f}")
    print(f"  总变化: {trend_full:.6e} ({trend_full/abs_full[-1]*100:.4f}%)")

    # 检查是否有系统性趋势
    if abs(trend_full) > 1e-4:
        print(f"\n  ⚠️  注意：完整版本显示系统性趋势（变化 > 0.01%）")
        print(f"  外推到 ε→0 可能需要更小的 epsilon 值")
    else:
        print(f"\n  ✅ 完整版本在测试范围内稳定")
