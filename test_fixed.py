#!/usr/bin/env python3
"""
测试修复后的径向流代码，绘制 σ 的演化图
"""

import numpy as np
import matplotlib.pyplot as plt
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

# 设置参数
bh = SchwarzschildParams(r0=2.0)
omega_test = 0.5
wave = WaveParams(omega=omega_test, l=2)

print(f"测试参数：ω = {omega_test}, l = 2")
print("=" * 60)

# 创建求解器并计算
solver = RadialFlowSolver(bh, wave)  
result = solver.solve(r_max=200.0)

print(f"积分成功：{result['success']}")
print(f"反射系数 A/B = {result['reflection_coeff']}")
print(f"|A/B| = {np.abs(result['reflection_coeff']):.6f}")
print()

# 提取数据
r = result['r']
sigma = result['sigma']
r0 = bh.r0

# 创建图像
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 图1：|sigma| vs r
axes[0, 0].plot(r/r0, np.abs(sigma), 'b-', linewidth=1.5)
axes[0, 0].set_xlabel(r'$r/r_0$', fontsize=12)
axes[0, 0].set_ylabel(r'$|\sigma|$', fontsize=12)
axes[0, 0].set_title(r'$|\sigma|$ vs $r$ (修复后)', fontsize=13)
axes[0, 0].grid(True, alpha=0.3)

# 图2：Re(sigma) 和 Im(sigma) vs r
axes[0, 1].plot(r/r0, sigma.real, 'r-', label=r'Re($\sigma$)', linewidth=1.5)
axes[0, 1].plot(r/r0, sigma.imag, 'b-', label=r'Im($\sigma$)', linewidth=1.5)
axes[0, 1].set_xlabel(r'$r/r_0$', fontsize=12)
axes[0, 1].set_ylabel(r'$\sigma$', fontsize=12)
axes[0, 1].set_title(r'$\sigma$ components vs $r$ (修复后)', fontsize=13)
axes[0, 1].legend(fontsize=11)
axes[0, 1].grid(True, alpha=0.3)

# 图3：复平面轨迹（全部点）
axes[1, 0].plot(sigma.real, sigma.imag, 'b-', alpha=0.6, linewidth=1)
axes[1, 0].scatter([sigma.real[0]], [sigma.imag[0]], c='g', s=100,
                   label='Start (horizon)', zorder=5, marker='o')
axes[1, 0].scatter([sigma.real[-1]], [sigma.imag[-1]], c='r', s=100,
                   label='End (infinity)', zorder=5, marker='s')
axes[1, 0].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[1, 0].axvline(x=0, color='k', linestyle='--', alpha=0.3)
axes[1, 0].set_xlabel(r'Re($\sigma$)', fontsize=12)
axes[1, 0].set_ylabel(r'Im($\sigma$)', fontsize=12)
axes[1, 0].set_title(r'$\sigma$ in complex plane (修复后)', fontsize=13)
axes[1, 0].legend(fontsize=10)
axes[1, 0].grid(True, alpha=0.3)
axes[1, 0].axis('equal')

# 图4：复平面轨迹（最后1000个点，放大）
n_last = min(1000, len(sigma))
sigma_last = sigma[-n_last:]
axes[1, 1].plot(sigma_last.real, sigma_last.imag, 'r-', linewidth=2)
axes[1, 1].scatter([sigma_last[0].real], [sigma_last[0].imag], c='g', s=100,
                   label='Start', zorder=5, marker='o')
axes[1, 1].scatter([sigma_last[-1].real], [sigma_last[-1].imag], c='b', s=100,
                   label='End', zorder=5, marker='s')
axes[1, 1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[1, 1].axvline(x=0, color='k', linestyle='--', alpha=0.3)
axes[1, 1].set_xlabel(r'Re($\sigma$)', fontsize=12)
axes[1, 1].set_ylabel(r'Im($\sigma$)', fontsize=12)
axes[1, 1].set_title(f'Zoomed: Last {n_last} points (修复后)', fontsize=13)
axes[1, 1].legend(fontsize=10)
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].axis('equal')

plt.tight_layout()

# 保存图像
save_path = f'sigma_fixed_omega{omega_test}.png'
plt.savefig(save_path, dpi=150, bbox_inches='tight')
print(f"图像已保存到 {save_path}")

# 显示图像
plt.show()
