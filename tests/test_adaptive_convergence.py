#!/usr/bin/env python3
"""
测试不同频率下的自适应收敛行为
"""

import numpy as np
import matplotlib.pyplot as plt
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

# 设置参数
bh = SchwarzschildParams(r0=2.0)

# 测试频率范围（从小到大）
omegas = np.array([0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0])

print("=" * 80)
print("自适应收敛测试：ω 从小到大")
print("=" * 80)
print()
print(f"{'ω':>6} | {'r_final/r0':>12} | {'n_points':>10} | {'收敛':>6} | {'|A/B|':>10} | {'A/B':>25}")
print("-" * 80)

results = []

for omega in omegas:
    wave = WaveParams(omega=omega, l=2)
    solver = RadialFlowSolver(bh, wave)

    # 使用自适应收敛
    result = solver.solve(
        r_max=200.0,
        convergence_threshold=0.05,  # 5% 阈值（放宽一点）
        n_periods=3,  # 检测3个周期
        max_r_factor=10.0  # 最多扩展到 10*r_max
    )

    r = result['r']
    r_final = r[-1] / bh.r0
    n_points = len(r)
    converged = result['success']
    A_over_B = result['reflection_coeff']
    abs_A_over_B = np.abs(A_over_B)

    # 保存结果
    results.append({
        'omega': omega,
        'r_final': r_final,
        'n_points': n_points,
        'converged': converged,
        'A_over_B': A_over_B,
        'abs_A_over_B': abs_A_over_B
    })

    # 打印结果
    converged_str = "是" if converged else "否"
    print(f"{omega:6.1f} | {r_final:12.1f} | {n_points:10d} | {converged_str:>6} | {abs_A_over_B:10.6f} | {A_over_B.real:+.4f}{A_over_B.imag:+.4f}j")

print("=" * 80)
print()

# 绘制结果
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 图1：收敛的 r_final vs ω
omegas_arr = np.array([r['omega'] for r in results])
r_finals = np.array([r['r_final'] for r in results])
converged_mask = np.array([r['converged'] for r in results])

axes[0, 0].semilogy(omegas_arr[converged_mask], r_finals[converged_mask], 'go-',
                    label='Converged', markersize=8, linewidth=2)
axes[0, 0].semilogy(omegas_arr[~converged_mask], r_finals[~converged_mask], 'rx',
                    label='Not converged', markersize=10, markeredgewidth=2)
axes[0, 0].set_xlabel(r'$\omega$', fontsize=13)
axes[0, 0].set_ylabel(r'$r_{\rm final}/r_0$', fontsize=13)
axes[0, 0].set_title('Convergence radius vs frequency', fontsize=14)
axes[0, 0].legend(fontsize=11)
axes[0, 0].grid(True, alpha=0.3)

# 图2：|A/B| vs ω
abs_A_over_B_arr = np.array([r['abs_A_over_B'] for r in results])
axes[0, 1].plot(omegas_arr, abs_A_over_B_arr, 'bo-', linewidth=2, markersize=8)
axes[0, 1].set_xlabel(r'$\omega$', fontsize=13)
axes[0, 1].set_ylabel(r'$|A/B|$', fontsize=13)
axes[0, 1].set_title('Reflection coefficient magnitude', fontsize=14)
axes[0, 1].grid(True, alpha=0.3)
axes[0, 1].set_ylim([0, 1.1])

# 图3：积分点数 vs ω
n_points_arr = np.array([r['n_points'] for r in results])
axes[1, 0].semilogy(omegas_arr, n_points_arr, 'mo-', linewidth=2, markersize=8)
axes[1, 0].set_xlabel(r'$\omega$', fontsize=13)
axes[1, 0].set_ylabel('Number of points', fontsize=13)
axes[1, 0].set_title('Integration points vs frequency', fontsize=14)
axes[1, 0].grid(True, alpha=0.3)

# 图4：A/B 在复平面
A_over_B_arr = np.array([r['A_over_B'] for r in results])
axes[1, 1].plot(A_over_B_arr.real, A_over_B_arr.imag, 'o-', linewidth=2, markersize=8)
for i, omega in enumerate(omegas_arr):
    axes[1, 1].annotate(f'{omega:.1f}',
                       (A_over_B_arr[i].real, A_over_B_arr[i].imag),
                       xytext=(5, 5), textcoords='offset points', fontsize=9)
axes[1, 1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[1, 1].axvline(x=0, color='k', linestyle='--', alpha=0.3)
axes[1, 1].set_xlabel(r'Re($A/B$)', fontsize=13)
axes[1, 1].set_ylabel(r'Im($A/B$)', fontsize=13)
axes[1, 1].set_title('Reflection coefficient in complex plane', fontsize=14)
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].axis('equal')

plt.tight_layout()
plt.savefig('adaptive_convergence_test.png', dpi=150, bbox_inches='tight')
print("图像已保存到 adaptive_convergence_test.png")
plt.show()

print()
print("物理验证：")
print(f"  低频 (ω={omegas[0]:.1f}): |A/B| = {results[0]['abs_A_over_B']:.6f} (应接近1，完全反射)")
print(f"  高频 (ω={omegas[-1]:.1f}): |A/B| = {results[-1]['abs_A_over_B']:.6f} (应接近0，完全透射)")
