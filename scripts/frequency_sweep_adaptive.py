#!/usr/bin/env python3
"""
频率扫描：测试自适应收敛
"""

import numpy as np
import matplotlib.pyplot as plt
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)

# 测试频率范围
omegas = np.array([0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0])

results = []

print("=" * 70)
print("频率扫描：自适应收敛测试")
print("=" * 70)
print()

for omega in omegas:
    wave = WaveParams(omega=omega, l=0)
    solver = RadialFlowSolver(bh, wave)

    result = solver.solve(
        r_max=200.0,
        convergence_threshold=0.1,  # 10% 阈值
        n_periods=5,  # 使用5个周期
        max_r_factor=10.0
    )

    r_final = result['r'][-1] / bh.r0
    n_points = len(result['r'])
    converged = result['success']
    reflection_coeff = result['reflection_coeff']
    abs_reflection = np.abs(reflection_coeff)

    results.append({
        'omega': omega,
        'r_final': r_final,
        'n_points': n_points,
        'converged': converged,
        'reflection_coeff': reflection_coeff,
        'abs_reflection': abs_reflection
    })

    print(f"ω = {omega:5.2f}: r/r0 = {r_final:7.1f}, "
          f"points = {n_points:6d}, converged = {converged}, "
          f"|A/B| = {abs_reflection:.6f}")

print()
print("=" * 70)
print("扫描完成")
print("=" * 70)

# 绘图
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 图1：|A/B| vs ω
axes[0, 0].plot([r['omega'] for r in results],
                [r['abs_reflection'] for r in results],
                'bo-', linewidth=2, markersize=8)
axes[0, 0].set_xlabel(r'$\omega$', fontsize=12)
axes[0, 0].set_ylabel(r'$|A/B|$', fontsize=12)
axes[0, 0].set_title('Reflection Coefficient vs Frequency', fontsize=13)
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_ylim([0, 1.1])

# 图2：收敛位置 r/r0 vs ω
axes[0, 1].plot([r['omega'] for r in results],
                [r['r_final'] for r in results],
                'ro-', linewidth=2, markersize=8)
axes[0, 1].set_xlabel(r'$\omega$', fontsize=12)
axes[0, 1].set_ylabel(r'$r_{\rm final}/r_0$', fontsize=12)
axes[0, 1].set_title('Convergence Position vs Frequency', fontsize=13)
axes[0, 1].grid(True, alpha=0.3)

# 图3：积分点数 vs ω
axes[1, 0].plot([r['omega'] for r in results],
                [r['n_points'] for r in results],
                'go-', linewidth=2, markersize=8)
axes[1, 0].set_xlabel(r'$\omega$', fontsize=12)
axes[1, 0].set_ylabel('Number of Points', fontsize=12)
axes[1, 0].set_title('Integration Points vs Frequency', fontsize=13)
axes[1, 0].grid(True, alpha=0.3)

# 图4：收敛状态
converged_omegas = [r['omega'] for r in results if r['converged']]
not_converged_omegas = [r['omega'] for r in results if not r['converged']]
converged_abs = [r['abs_reflection'] for r in results if r['converged']]
not_converged_abs = [r['abs_reflection'] for r in results if not r['converged']]

axes[1, 1].plot(converged_omegas, converged_abs,
                'go', markersize=10, label='Converged')
if not_converged_omegas:
    axes[1, 1].plot(not_converged_omegas, not_converged_abs,
                    'rx', markersize=10, label='Not Converged')
axes[1, 1].set_xlabel(r'$\omega$', fontsize=12)
axes[1, 1].set_ylabel(r'$|A/B|$', fontsize=12)
axes[1, 1].set_title('Convergence Status', fontsize=13)
axes[1, 1].legend(fontsize=11)
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].set_ylim([0, 1.1])

plt.tight_layout()
plt.savefig('frequency_sweep_adaptive.png', dpi=150, bbox_inches='tight')
print("\n图像已保存到 frequency_sweep_adaptive.png")
# plt.show()  # Commented out to avoid blocking
