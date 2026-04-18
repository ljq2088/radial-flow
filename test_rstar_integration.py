#!/usr/bin/env python3
"""
测试在r*坐标积分并推进到r=2000的结果
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 添加python目录到路径
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

def test_rstar_integration():
    """测试r*坐标积分"""

    bh = SchwarzschildParams(r0=2.0)
    omega = 0.2
    l = 0

    print("=" * 70)
    print("测试在r*坐标积分并推进到r=2000")
    print("=" * 70)
    print(f"参数: ω = {omega}, l = {l}")
    print()

    # 测试简化版本
    print("-" * 70)
    print("简化版本 (V = l(l+1)/r²)")
    print("-" * 70)

    # 在r坐标积分
    wave_simple_r = WaveParams(omega=omega, l=l, use_subleading_term=False, integrate_in_rstar=False)
    solver_simple_r = RadialFlowSolver(bh, wave_simple_r)
    result_simple_r = solver_simple_r.solve(r_max=2000.0, epsilon=1e-6, n_points=50000)

    print(f"在r坐标积分:")
    print(f"  最终r/r₀ = {result_simple_r['r'][-1]/bh.r0:.2f}")
    print(f"  sigma_final = {result_simple_r['sigma'][-1]}")
    print(f"  A/B = {result_simple_r['reflection_coeff']}")
    print(f"  |A/B| = {abs(result_simple_r['reflection_coeff']):.6f}")
    print(f"  arg(A/B) = {np.angle(result_simple_r['reflection_coeff'])/np.pi:.6f}π")
    print(f"  converged = {result_simple_r['success']}")
    print()

    # 在r*坐标积分
    wave_simple_rstar = WaveParams(omega=omega, l=l, use_subleading_term=False, integrate_in_rstar=True)
    solver_simple_rstar = RadialFlowSolver(bh, wave_simple_rstar)
    result_simple_rstar = solver_simple_rstar.solve(r_max=2000.0, epsilon=1e-6, n_points=50000)

    print(f"在r*坐标积分:")
    print(f"  最终r/r₀ = {result_simple_rstar['r'][-1]/bh.r0:.2f}")
    print(f"  sigma_final = {result_simple_rstar['sigma'][-1]}")
    print(f"  A/B = {result_simple_rstar['reflection_coeff']}")
    print(f"  |A/B| = {abs(result_simple_rstar['reflection_coeff']):.6f}")
    print(f"  arg(A/B) = {np.angle(result_simple_rstar['reflection_coeff'])/np.pi:.6f}π")
    print(f"  converged = {result_simple_rstar['success']}")
    print()

    # 比较
    diff_mag = abs(result_simple_rstar['reflection_coeff']) - abs(result_simple_r['reflection_coeff'])
    diff_phase = (np.angle(result_simple_rstar['reflection_coeff']) - np.angle(result_simple_r['reflection_coeff'])) / np.pi
    print(f"差异:")
    print(f"  Δ|A/B| = {diff_mag:.6f} ({100*diff_mag/abs(result_simple_r['reflection_coeff']):.3f}%)")
    print(f"  Δarg(A/B) = {diff_phase:.6f}π")
    print()

    # 测试完整版本
    print("-" * 70)
    print("完整版本 (V = l(l+1)/r² + f'/r)")
    print("-" * 70)

    # 在r坐标积分
    wave_full_r = WaveParams(omega=omega, l=l, use_subleading_term=True, integrate_in_rstar=False)
    solver_full_r = RadialFlowSolver(bh, wave_full_r)
    result_full_r = solver_full_r.solve(r_max=2000.0, epsilon=1e-6, n_points=50000)

    print(f"在r坐标积分:")
    print(f"  最终r/r₀ = {result_full_r['r'][-1]/bh.r0:.2f}")
    print(f"  sigma_final = {result_full_r['sigma'][-1]}")
    print(f"  A/B = {result_full_r['reflection_coeff']}")
    print(f"  |A/B| = {abs(result_full_r['reflection_coeff']):.6f}")
    print(f"  arg(A/B) = {np.angle(result_full_r['reflection_coeff'])/np.pi:.6f}π")
    print(f"  converged = {result_full_r['success']}")
    print()

    # 在r*坐标积分
    wave_full_rstar = WaveParams(omega=omega, l=l, use_subleading_term=True, integrate_in_rstar=True)
    solver_full_rstar = RadialFlowSolver(bh, wave_full_rstar)
    result_full_rstar = solver_full_rstar.solve(r_max=2000.0, epsilon=1e-6, n_points=50000)

    print(f"在r*坐标积分:")
    print(f"  最终r/r₀ = {result_full_rstar['r'][-1]/bh.r0:.2f}")
    print(f"  sigma_final = {result_full_rstar['sigma'][-1]}")
    print(f"  A/B = {result_full_rstar['reflection_coeff']}")
    print(f"  |A/B| = {abs(result_full_rstar['reflection_coeff']):.6f}")
    print(f"  arg(A/B) = {np.angle(result_full_rstar['reflection_coeff'])/np.pi:.6f}π")
    print(f"  converged = {result_full_rstar['success']}")
    print()

    # 比较
    diff_mag = abs(result_full_rstar['reflection_coeff']) - abs(result_full_r['reflection_coeff'])
    diff_phase = (np.angle(result_full_rstar['reflection_coeff']) - np.angle(result_full_r['reflection_coeff'])) / np.pi
    print(f"差异:")
    print(f"  Δ|A/B| = {diff_mag:.6f} ({100*diff_mag/abs(result_full_r['reflection_coeff']):.3f}%)")
    print(f"  Δarg(A/B) = {diff_phase:.6f}π")
    print()

    # 可视化
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 简化版本 - |σ|
    ax = axes[0, 0]
    r_simple_r = np.array(result_simple_r['r'])
    sigma_simple_r = np.array(result_simple_r['sigma'])
    r_simple_rstar = np.array(result_simple_rstar['r'])
    sigma_simple_rstar = np.array(result_simple_rstar['sigma'])

    ax.plot(r_simple_r/bh.r0, np.abs(sigma_simple_r), 'b-', label='r coordinate', linewidth=1)
    ax.plot(r_simple_rstar/bh.r0, np.abs(sigma_simple_rstar), 'r--', label='r* coordinate', linewidth=1)
    ax.set_xlabel('r/r₀')
    ax.set_ylabel('|σ|')
    ax.set_title('Simplified version: |σ| vs r')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 简化版本 - arg(σ)
    ax = axes[0, 1]
    ax.plot(r_simple_r/bh.r0, np.angle(sigma_simple_r)/np.pi, 'b-', label='r coordinate', linewidth=1)
    ax.plot(r_simple_rstar/bh.r0, np.angle(sigma_simple_rstar)/np.pi, 'r--', label='r* coordinate', linewidth=1)
    ax.set_xlabel('r/r₀')
    ax.set_ylabel('arg(σ)/π')
    ax.set_title('Simplified version: arg(σ) vs r')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 完整版本 - |σ|
    ax = axes[1, 0]
    r_full_r = np.array(result_full_r['r'])
    sigma_full_r = np.array(result_full_r['sigma'])
    r_full_rstar = np.array(result_full_rstar['r'])
    sigma_full_rstar = np.array(result_full_rstar['sigma'])

    ax.plot(r_full_r/bh.r0, np.abs(sigma_full_r), 'b-', label='r coordinate', linewidth=1)
    ax.plot(r_full_rstar/bh.r0, np.abs(sigma_full_rstar), 'r--', label='r* coordinate', linewidth=1)
    ax.set_xlabel('r/r₀')
    ax.set_ylabel('|σ|')
    ax.set_title('Full version: |σ| vs r')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 完整版本 - arg(σ)
    ax = axes[1, 1]
    ax.plot(r_full_r/bh.r0, np.angle(sigma_full_r)/np.pi, 'b-', label='r coordinate', linewidth=1)
    ax.plot(r_full_rstar/bh.r0, np.angle(sigma_full_rstar)/np.pi, 'r--', label='r* coordinate', linewidth=1)
    ax.set_xlabel('r/r₀')
    ax.set_ylabel('arg(σ)/π')
    ax.set_title('Full version: arg(σ) vs r')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('figures/rstar_integration_comparison.png', dpi=150, bbox_inches='tight')
    print(f"图像已保存到: figures/rstar_integration_comparison.png")

if __name__ == "__main__":
    test_rstar_integration()
