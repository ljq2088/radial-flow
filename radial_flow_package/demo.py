#!/usr/bin/env python3
"""
Radial Flow 计算演示脚本
======================

计算 Schwarzschild 黑洞背景下标量场的反射系数 A/B
对比简化版本和完整版本的结果，并生成可视化图表
"""

import sys
import os

# 添加当前目录到 Python 路径（使用相对路径）
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver
import numpy as np
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt
from datetime import datetime

# ============================================================================
# 参数设置
# ============================================================================

# 黑洞参数
bh = SchwarzschildParams(r0=2.0)  # Schwarzschild 半径 (M=1)

# 波参数
omega = 0.2  # 频率
l = 0        # 角动量量子数

# 求解参数
r_max = 200.0
n_points = 10000
convergence_threshold = 0.01
n_periods = 2
max_r_factor = 2.0

# ============================================================================
# 计算简化版本和完整版本
# ============================================================================

print("=" * 80)
print("Schwarzschild 黑洞标量场散射计算")
print("=" * 80)
print(f"计算时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"\n黑洞参数:")
print(f"  Schwarzschild 半径: r₀ = {bh.r0}")
print(f"  质量: M = {bh.r0/2}")
print(f"\n波参数:")
print(f"  频率: ω = {omega}")
print(f"  角动量量子数: l = {l}")
print(f"\n求解参数:")
print(f"  最小积分终点: r_max = {r_max}")
print(f"  输出点数: n_points = {n_points}")
print(f"  收敛阈值: {convergence_threshold*100}%")
print(f"  检测周期数: {n_periods}")
print(f"  最大扩展倍数: {max_r_factor}")
print("\n" + "-" * 80)

# 简化版本: V = l(l+1)/r²
print("\n[1/2] 计算简化版本 (V = l(l+1)/r²)...")
wave_simple = WaveParams(omega=omega, l=l, use_subleading_term=False)
solver_simple = RadialFlowSolver(bh, wave_simple)
result_simple = solver_simple.solve(
    r_max=r_max,
    n_points=n_points,
    convergence_threshold=convergence_threshold,
    n_periods=n_periods,
    max_r_factor=max_r_factor
)

# 完整版本: V = l(l+1)/r² + f'/r
print("[2/2] 计算完整版本 (V = l(l+1)/r² + f'/r)...")
wave_full = WaveParams(omega=omega, l=l, use_subleading_term=True)
solver_full = RadialFlowSolver(bh, wave_full)
result_full = solver_full.solve(
    r_max=r_max,
    n_points=n_points,
    convergence_threshold=convergence_threshold,
    n_periods=n_periods,
    max_r_factor=max_r_factor
)

# ============================================================================
# 文本输出结果
# ============================================================================

print("\n" + "=" * 80)
print("计算结果")
print("=" * 80)

# 简化版本结果
A_B_simple = result_simple['reflection_coeff']
print(f"\n【简化版本】V = l(l+1)/r²")
print(f"  势函数: 仅包含角动量势垒")
print(f"  收敛状态: {result_simple['success']}")
print(f"  积分终点: r = {result_simple['r'][-1]:.2f}")
print(f"  反射系数: A/B = {A_B_simple.real:.6f} + {A_B_simple.imag:.6f}i")
print(f"  模: |A/B| = {abs(A_B_simple):.6f}")
print(f"  相位: arg(A/B) = {np.angle(A_B_simple):.6f} rad")

# 完整版本结果
A_B_full = result_full['reflection_coeff']
print(f"\n【完整版本】V = l(l+1)/r² + f'/r")
print(f"  势函数: 包含角动量势垒 + 时空曲率修正")
print(f"  收敛状态: {result_full['success']}")
print(f"  积分终点: r = {result_full['r'][-1]:.2f}")
print(f"  反射系数: A/B = {A_B_full.real:.6f} + {A_B_full.imag:.6f}i")
print(f"  模: |A/B| = {abs(A_B_full):.6f}")
print(f"  相位: arg(A/B) = {np.angle(A_B_full):.6f} rad")

# 对比分析
diff_abs = abs(abs(A_B_full) - abs(A_B_simple))
diff_percent = diff_abs / abs(A_B_simple) * 100
print(f"\n【对比分析】")
print(f"  |A/B| 差异: {diff_abs:.6f}")
print(f"  相对差异: {diff_percent:.2f}%")
print(f"  物理解释: 次领头项 f'/r 对反射系数有 {diff_percent:.1f}% 的修正")

# ============================================================================
# 可视化
# ============================================================================

print("\n" + "-" * 80)
print("生成可视化图表...")

fig = plt.figure(figsize=(16, 10))

# 子图1: |sigma| vs r (简化版本)
ax1 = plt.subplot(2, 3, 1)
r_simple = np.array(result_simple['r'])
sigma_simple = np.array(result_simple['sigma'])
ax1.plot(r_simple/bh.r0, np.abs(sigma_simple), 'b-', linewidth=1.5)
ax1.set_xlabel(r'$r/r_0$', fontsize=12)
ax1.set_ylabel(r'$|\sigma|$', fontsize=12)
ax1.set_title(r'Simplified: $|\sigma|$ vs $r$', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(left=1)

# 子图2: |sigma| vs r (完整版本)
ax2 = plt.subplot(2, 3, 2)
r_full = np.array(result_full['r'])
sigma_full = np.array(result_full['sigma'])
ax2.plot(r_full/bh.r0, np.abs(sigma_full), 'r-', linewidth=1.5)
ax2.set_xlabel(r'$r/r_0$', fontsize=12)
ax2.set_ylabel(r'$|\sigma|$', fontsize=12)
ax2.set_title(r'Full: $|\sigma|$ vs $r$', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(left=1)

# 子图3: 对比 |sigma|
ax3 = plt.subplot(2, 3, 3)
ax3.plot(r_simple/bh.r0, np.abs(sigma_simple), 'b-', linewidth=1.5, label='Simplified', alpha=0.7)
ax3.plot(r_full/bh.r0, np.abs(sigma_full), 'r-', linewidth=1.5, label='Full', alpha=0.7)
ax3.set_xlabel(r'$r/r_0$', fontsize=12)
ax3.set_ylabel(r'$|\sigma|$', fontsize=12)
ax3.set_title(r'Comparison: $|\sigma|$ vs $r$', fontsize=13, fontweight='bold')
ax3.legend(fontsize=11)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(left=1)

# 子图4: 复平面轨迹 (简化版本)
ax4 = plt.subplot(2, 3, 4)
ax4.plot(sigma_simple.real, sigma_simple.imag, 'b-', linewidth=1, alpha=0.6)
ax4.scatter([sigma_simple.real[0]], [sigma_simple.imag[0]], c='green', s=100,
            label='Start (horizon)', zorder=5, marker='o')
ax4.scatter([sigma_simple.real[-1]], [sigma_simple.imag[-1]], c='red', s=100,
            label='End (infinity)', zorder=5, marker='s')
ax4.set_xlabel(r'Re($\sigma$)', fontsize=12)
ax4.set_ylabel(r'Im($\sigma$)', fontsize=12)
ax4.set_title('Simplified: Complex plane', fontsize=13, fontweight='bold')
ax4.legend(fontsize=10)
ax4.grid(True, alpha=0.3)
ax4.axis('equal')

# 子图5: 复平面轨迹 (完整版本)
ax5 = plt.subplot(2, 3, 5)
ax5.plot(sigma_full.real, sigma_full.imag, 'r-', linewidth=1, alpha=0.6)
ax5.scatter([sigma_full.real[0]], [sigma_full.imag[0]], c='green', s=100,
            label='Start (horizon)', zorder=5, marker='o')
ax5.scatter([sigma_full.real[-1]], [sigma_full.imag[-1]], c='red', s=100,
            label='End (infinity)', zorder=5, marker='s')
ax5.set_xlabel(r'Re($\sigma$)', fontsize=12)
ax5.set_ylabel(r'Im($\sigma$)', fontsize=12)
ax5.set_title('Full: Complex plane', fontsize=13, fontweight='bold')
ax5.legend(fontsize=10)
ax5.grid(True, alpha=0.3)
ax5.axis('equal')

# 子图6: 结果总结
ax6 = plt.subplot(2, 3, 6)
ax6.axis('off')

summary_text = f"""
Parameters
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Black hole: r0 = {bh.r0}, M = {bh.r0/2}
Wave: omega = {omega}, l = {l}

Results
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Simplified (V = l(l+1)/r²)
  |A/B| = {abs(A_B_simple):.6f}
  arg(A/B) = {np.angle(A_B_simple):.4f} rad

Full (V = l(l+1)/r² + f'/r)
  |A/B| = {abs(A_B_full):.6f}
  arg(A/B) = {np.angle(A_B_full):.4f} rad

Analysis
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Absolute diff: {diff_abs:.6f}
Relative diff: {diff_percent:.2f}%

Physical meaning
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Subleading term f'/r causes
~{diff_percent:.1f}% correction to the
reflection coefficient
"""

ax6.text(0.1, 0.5, summary_text, transform=ax6.transAxes,
         fontsize=11, verticalalignment='center',
         family='monospace', bbox=dict(boxstyle='round',
         facecolor='wheat', alpha=0.3))

plt.suptitle(f'Schwarzschild Black Hole Scalar Scattering (omega={omega}, l={l})',
             fontsize=15, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])

# 保存图表
output_dir = os.path.join(current_dir, 'figures')
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'radial_flow_comparison.png')
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"图表已保存: {output_file}")

# 不显示图表窗口（适用于无头环境）
# plt.show()

print("\n" + "=" * 80)
print("计算完成！")
print("=" * 80)
