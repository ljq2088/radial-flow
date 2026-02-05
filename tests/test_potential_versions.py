#!/usr/bin/env python3
"""
测试两种势函数版本对 l=0, ω=0.2 的影响

比较：
1. 完整版本：V = l(l+1)/r² + f'/r  (.tex 第52-86行)
2. 简化版本：V = l(l+1)/r²        (.tex 第87-108行)
"""

import sys
sys.path.insert(0, '/home/ljq/code/radial_flow')

from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver
import matplotlib.pyplot as plt
import numpy as np

# 黑洞参数
bh = SchwarzschildParams(r0=2.0)

# 测试参数
omega = 0.2
l = 0

print("=" * 60)
print(f"测试两种势函数版本：l={l}, ω={omega}")
print("=" * 60)

# ===== 测试简化版本 =====
print("\n[1] 简化版本: V = l(l+1)/r²")
print("-" * 60)

wave_simple = WaveParams(omega=omega, l=l, use_full_potential=False)
solver_simple = RadialFlowSolver(bh, wave_simple)

result_simple = solver_simple.solve(
    r_max=200.0,
    epsilon=1e-6,
    n_points=100000,
    convergence_threshold=0.05,  # 放宽阈值，简化版本收敛性较差
    n_periods=2,
    max_r_factor=2.0  # 限制最大扩展倍数
)

print(f"收敛状态: {result_simple['success']}")
print(f"sigma_final = {result_simple['sigma'][-1]:.6f}")
print(f"反射系数 |A/B| = {abs(result_simple['reflection_coeff']):.6f}")
print(f"相位 arg(A/B) = {np.angle(result_simple['reflection_coeff']):.6f} rad")
print(f"实际 r_max = {result_simple['r'][-1]:.2f}")

# ===== 测试完整版本 =====
print("\n[2] 完整版本: V = l(l+1)/r² + f'/r")
print("-" * 60)

wave_full = WaveParams(omega=omega, l=l, use_full_potential=True)
solver_full = RadialFlowSolver(bh, wave_full)

result_full = solver_full.solve(
    r_max=200.0,
    epsilon=1e-6,
    n_points=100000,
    convergence_threshold=0.01,
    n_periods=2,
    max_r_factor=2.0  # 限制最大扩展倍数
)

print(f"收敛状态: {result_full['success']}")
print(f"sigma_final = {result_full['sigma'][-1]:.6f}")
print(f"反射系数 |A/B| = {abs(result_full['reflection_coeff']):.6f}")
print(f"相位 arg(A/B) = {np.angle(result_full['reflection_coeff']):.6f} rad")
print(f"实际 r_max = {result_full['r'][-1]:.2f}")

# ===== 对比分析 =====
print("\n[3] 对比分析")
print("=" * 60)

diff_amplitude = abs(abs(result_full['reflection_coeff']) - abs(result_simple['reflection_coeff']))
diff_phase = abs(np.angle(result_full['reflection_coeff']) - np.angle(result_simple['reflection_coeff']))

print(f"振幅差异: Δ|A/B| = {diff_amplitude:.6f}")
print(f"相位差异: Δarg(A/B) = {diff_phase:.6f} rad")

# ===== 可视化 =====
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 提取数据
r_simple = np.array(result_simple['r'])
sigma_simple = np.array(result_simple['sigma'])
r_full = np.array(result_full['r'])
sigma_full = np.array(result_full['sigma'])

# 绘制 sigma 实部
ax = axes[0, 0]
ax.plot(r_simple, np.real(sigma_simple), 'b-', label='简化版本', alpha=0.7)
ax.plot(r_full, np.real(sigma_full), 'r--', label='完整版本', alpha=0.7)
ax.axhline(y=0, color='k', linestyle=':', alpha=0.3)
ax.set_xlabel('r / r₀')
ax.set_ylabel('Re(σ)')
ax.set_title('σ 的实部')
ax.legend()
ax.grid(True, alpha=0.3)

# 绘制 sigma 虚部
ax = axes[0, 1]
ax.plot(r_simple, np.imag(sigma_simple), 'b-', label='简化版本', alpha=0.7)
ax.plot(r_full, np.imag(sigma_full), 'r--', label='完整版本', alpha=0.7)
ax.axhline(y=0, color='k', linestyle=':', alpha=0.3)
ax.set_xlabel('r / r₀')
ax.set_ylabel('Im(σ)')
ax.set_title('σ 的虚部')
ax.legend()
ax.grid(True, alpha=0.3)

# 复平面轨迹
ax = axes[1, 0]
ax.plot(np.real(sigma_simple), np.imag(sigma_simple), 'b-', label='简化版本', alpha=0.7, linewidth=1)
ax.plot(np.real(sigma_full), np.imag(sigma_full), 'r--', label='完整版本', alpha=0.7, linewidth=1)
ax.plot([0], [0], 'go', markersize=8, label='起点 (σ₀=0)')
ax.plot([np.real(sigma_simple[-1])], [np.imag(sigma_simple[-1])], 'bs', markersize=8, label='简化终点')
ax.plot([np.real(sigma_full[-1])], [np.imag(sigma_full[-1])], 'r^', markersize=8, label='完整终点')
ax.set_xlabel('Re(σ)')
ax.set_ylabel('Im(σ)')
ax.set_title('复平面轨迹')
ax.legend(loc='best', fontsize=8)
ax.grid(True, alpha=0.3)
ax.axis('equal')

# |sigma| 随 r 变化
ax = axes[1, 1]
ax.plot(r_simple, np.abs(sigma_simple), 'b-', label='简化版本', alpha=0.7)
ax.plot(r_full, np.abs(sigma_full), 'r--', label='完整版本', alpha=0.7)
ax.set_xlabel('r / r₀')
ax.set_ylabel('|σ|')
ax.set_title('|σ| 随半径变化')
ax.legend()
ax.grid(True, alpha=0.3)

plt.suptitle(f'两种势函数版本对比: l={l}, ω={omega}', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/home/ljq/code/radial_flow/figures/potential_versions_comparison.png', dpi=150, bbox_inches='tight')
print(f"\n图片已保存至: figures/potential_versions_comparison.png")
plt.show()
