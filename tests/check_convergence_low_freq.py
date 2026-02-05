#!/usr/bin/env python3
"""
检查低频情况的收敛行为
"""

import numpy as np
import matplotlib.pyplot as plt
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

# 测试 ω=0.1
bh = SchwarzschildParams(r0=2.0)
omega = 0.1
wave = WaveParams(omega=omega, l=2)
solver = RadialFlowSolver(bh, wave)

print(f"测试 ω = {omega}")
print("=" * 60)

# 运行求解
result = solver.solve(
    r_max=200.0,
    convergence_threshold=0.1,  # 10% 阈值
    n_periods=2,
    max_r_factor=10.0
)

r = result['r']
sigma = result['sigma']
r0 = bh.r0

print(f"最终 r/r0 = {r[-1]/r0:.1f}")
print(f"积分点数 = {len(r)}")
print(f"收敛状态 = {result['success']}")
print(f"反射系数 A/B = {result['reflection_coeff']}")
print(f"|A/B| = {np.abs(result['reflection_coeff']):.6f}")
print()

# 分析最后一段的震荡
n_check = min(5000, len(sigma))
sigma_last = sigma[-n_check:]
r_last = r[-n_check:]

# 计算平均值
sigma_mean = np.mean(sigma_last)
print(f"最后 {n_check} 点的平均值: {sigma_mean}")
print(f"|平均值| = {np.abs(sigma_mean):.6f}")
print()

# 计算震荡幅度
deviations = np.abs(sigma_last - sigma_mean)
max_deviation = np.max(deviations)
mean_deviation = np.mean(deviations)
relative_oscillation = max_deviation / np.abs(sigma_mean) if np.abs(sigma_mean) > 0 else np.inf

print(f"最大偏差 = {max_deviation:.6f}")
print(f"平均偏差 = {mean_deviation:.6f}")
print(f"相对震荡幅度 = {relative_oscillation:.6f} ({relative_oscillation*100:.2f}%)")
print()

# 检查震荡周期
period_estimate = 2 * np.pi / omega
dr = np.mean(np.diff(r_last))
points_per_period = period_estimate / dr
print(f"估计震荡周期 = {period_estimate:.2f}")
print(f"每周期点数 ≈ {points_per_period:.0f}")
print()

# 绘图
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# 图1：|σ| vs r（全部）
axes[0, 0].plot(r/r0, np.abs(sigma), 'b-', linewidth=1)
axes[0, 0].set_xlabel(r'$r/r_0$')
axes[0, 0].set_ylabel(r'$|\sigma|$')
axes[0, 0].set_title(f'$|\sigma|$ vs $r$ (ω={omega})')
axes[0, 0].grid(True, alpha=0.3)

# 图2：|σ| vs r（最后部分）
axes[0, 1].plot(r_last/r0, np.abs(sigma_last), 'b-', linewidth=1)
axes[0, 1].axhline(y=np.abs(sigma_mean), color='r', linestyle='--', label='Mean')
axes[0, 1].set_xlabel(r'$r/r_0$')
axes[0, 1].set_ylabel(r'$|\sigma|$')
axes[0, 1].set_title(f'$|\sigma|$ vs $r$ (last {n_check} points)')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# 图3：Re(σ) 和 Im(σ) vs r（最后部分）
axes[0, 2].plot(r_last/r0, sigma_last.real, 'r-', label='Re', linewidth=1)
axes[0, 2].plot(r_last/r0, sigma_last.imag, 'b-', label='Im', linewidth=1)
axes[0, 2].axhline(y=sigma_mean.real, color='r', linestyle='--', alpha=0.5)
axes[0, 2].axhline(y=sigma_mean.imag, color='b', linestyle='--', alpha=0.5)
axes[0, 2].set_xlabel(r'$r/r_0$')
axes[0, 2].set_ylabel(r'$\sigma$')
axes[0, 2].set_title('Components (last part)')
axes[0, 2].legend()
axes[0, 2].grid(True, alpha=0.3)

# 图4：σ 在复平面（全部）
axes[1, 0].plot(sigma.real, sigma.imag, 'b-', alpha=0.6, linewidth=1)
axes[1, 0].scatter([sigma[0].real], [sigma[0].imag], c='g', s=100, marker='o', label='Start')
axes[1, 0].scatter([sigma[-1].real], [sigma[-1].imag], c='r', s=100, marker='s', label='End')
axes[1, 0].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[1, 0].axvline(x=0, color='k', linestyle='--', alpha=0.3)
axes[1, 0].set_xlabel(r'Re($\sigma$)')
axes[1, 0].set_ylabel(r'Im($\sigma$)')
axes[1, 0].set_title('Complex plane (all)')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)
axes[1, 0].axis('equal')

# 图5：σ 在复平面（最后部分）
axes[1, 1].plot(sigma_last.real, sigma_last.imag, 'r-', linewidth=2)
axes[1, 1].scatter([sigma_mean.real], [sigma_mean.imag], c='k', s=100, marker='x', label='Mean', zorder=5)
axes[1, 1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[1, 1].axvline(x=0, color='k', linestyle='--', alpha=0.3)
axes[1, 1].set_xlabel(r'Re($\sigma$)')
axes[1, 1].set_ylabel(r'Im($\sigma$)')
axes[1, 1].set_title(f'Complex plane (last {n_check} points)')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].axis('equal')

# 图6：偏差 vs r
axes[1, 2].plot(r_last/r0, deviations, 'g-', linewidth=1)
axes[1, 2].axhline(y=max_deviation, color='r', linestyle='--', label=f'Max = {max_deviation:.4f}')
axes[1, 2].axhline(y=mean_deviation, color='b', linestyle='--', label=f'Mean = {mean_deviation:.4f}')
axes[1, 2].set_xlabel(r'$r/r_0$')
axes[1, 2].set_ylabel(r'$|\sigma - \bar{\sigma}|$')
axes[1, 2].set_title('Deviation from mean')
axes[1, 2].legend()
axes[1, 2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'convergence_check_omega{omega}.png', dpi=150, bbox_inches='tight')
print(f"图像已保存到 convergence_check_omega{omega}.png")
plt.show()
