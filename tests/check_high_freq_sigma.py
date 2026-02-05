#!/usr/bin/env python3
"""
检查高频时 σ 的数值行为
"""

import numpy as np
import matplotlib.pyplot as plt
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)
omega = 3.0
wave = WaveParams(omega=omega, l=2)
solver = RadialFlowSolver(bh, wave)

print(f"测试 ω = {omega}")
print("=" * 60)

# 使用固定的 r_max，不使用收敛判据
result = solver.solve(
    r_max=200.0,
    convergence_threshold=1e10,  # 设置很大，强制积分到 r_max
    n_periods=2,
    max_r_factor=1.0  # 不扩展
)

r = result['r']
sigma = result['sigma']
r0 = bh.r0

print(f"积分点数 = {len(r)}")
print(f"最终 r/r0 = {r[-1]/r0:.1f}")
print()

# 检查 σ 的值
print("检查 σ 的数值：")
print(f"  σ[0] = {sigma[0]}")
print(f"  σ[100] = {sigma[100]}")
print(f"  σ[1000] = {sigma[1000]}")
print(f"  σ[-1000] = {sigma[-1000]}")
print(f"  σ[-100] = {sigma[-100]}")
print(f"  σ[-1] = {sigma[-1]}")
print()

# 检查是否有 nan 或 inf
has_nan = np.any(np.isnan(sigma))
has_inf = np.any(np.isinf(sigma))
print(f"包含 nan: {has_nan}")
print(f"包含 inf: {has_inf}")
print()

# 检查 |σ| 的范围
abs_sigma = np.abs(sigma)
print(f"|σ| 的范围：")
print(f"  最小值 = {np.min(abs_sigma):.6e}")
print(f"  最大值 = {np.max(abs_sigma):.6e}")
print(f"  平均值 = {np.mean(abs_sigma):.6e}")
print()

# 检查反射系数
print(f"反射系数 A/B = {result['reflection_coeff']}")
print(f"|A/B| = {np.abs(result['reflection_coeff'])}")
print()

# 绘图
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 图1：|σ| vs r
axes[0, 0].plot(r/r0, abs_sigma, 'b-', linewidth=1)
axes[0, 0].set_xlabel(r'$r/r_0$')
axes[0, 0].set_ylabel(r'$|\sigma|$')
axes[0, 0].set_title(f'$|\sigma|$ vs $r$ (ω={omega})')
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_yscale('log')

# 图2：Re(σ) 和 Im(σ) vs r
axes[0, 1].plot(r/r0, sigma.real, 'r-', label='Re', linewidth=1)
axes[0, 1].plot(r/r0, sigma.imag, 'b-', label='Im', linewidth=1)
axes[0, 1].set_xlabel(r'$r/r_0$')
axes[0, 1].set_ylabel(r'$\sigma$')
axes[0, 1].set_title('Components')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# 图3：σ 在复平面
axes[1, 0].plot(sigma.real, sigma.imag, 'b-', alpha=0.6, linewidth=1)
axes[1, 0].scatter([sigma[0].real], [sigma[0].imag], c='g', s=100, marker='o', label='Start')
axes[1, 0].scatter([sigma[-1].real], [sigma[-1].imag], c='r', s=100, marker='s', label='End')
axes[1, 0].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[1, 0].axvline(x=0, color='k', linestyle='--', alpha=0.3)
axes[1, 0].set_xlabel(r'Re($\sigma$)')
axes[1, 0].set_ylabel(r'Im($\sigma$)')
axes[1, 0].set_title('Complex plane')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)
axes[1, 0].axis('equal')

# 图4：最后一段的 |σ|
n_last = min(2000, len(sigma))
axes[1, 1].plot(r[-n_last:]/r0, abs_sigma[-n_last:], 'g-', linewidth=1)
axes[1, 1].set_xlabel(r'$r/r_0$')
axes[1, 1].set_ylabel(r'$|\sigma|$')
axes[1, 1].set_title(f'Last {n_last} points')
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].set_yscale('log')

plt.tight_layout()
plt.savefig(f'high_freq_omega{omega}_check.png', dpi=150, bbox_inches='tight')
print(f"图像已保存到 high_freq_omega{omega}_check.png")
