#!/usr/bin/env python3
"""
可视化 epsilon 外推行为，并测试极小的 epsilon 值
"""

import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 测试参数
bh = SchwarzschildParams(r0=2.0)
omega = 0.2
l = 0

# 扩展的 epsilon 值范围（包括极小值）
epsilon_values = [1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5,
                  5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9]

print("计算扩展的 epsilon 范围...")

# 测试简化版本
wave_simple = WaveParams(omega=omega, l=l, use_subleading_term=False)
solver_simple = RadialFlowSolver(bh, wave_simple)

results_simple = []
for i, eps in enumerate(epsilon_values):
    print(f"简化版本: {i+1}/{len(epsilon_values)} (ε={eps:.1e})", end='\r')
    try:
        result = solver_simple.solve(r_max=200.0, epsilon=eps, n_points=10000)
        A_B = result['reflection_coeff']
        results_simple.append({'epsilon': eps, 'abs': abs(A_B)})
    except:
        results_simple.append({'epsilon': eps, 'abs': None})

print()

# 测试完整版本
wave_full = WaveParams(omega=omega, l=l, use_subleading_term=True)
solver_full = RadialFlowSolver(bh, wave_full)

results_full = []
for i, eps in enumerate(epsilon_values):
    print(f"完整版本: {i+1}/{len(epsilon_values)} (ε={eps:.1e})", end='\r')
    try:
        result = solver_full.solve(r_max=200.0, epsilon=eps, n_points=10000)
        A_B = result['reflection_coeff']
        results_full.append({'epsilon': eps, 'abs': abs(A_B)})
    except:
        results_full.append({'epsilon': eps, 'abs': None})

print()
print("生成可视化...")

# 提取有效数据
eps_simple = [r['epsilon'] for r in results_simple if r['abs'] is not None]
abs_simple = [r['abs'] for r in results_simple if r['abs'] is not None]

eps_full = [r['epsilon'] for r in results_full if r['abs'] is not None]
abs_full = [r['abs'] for r in results_full if r['abs'] is not None]

# 创建图表
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 子图1: |A/B| vs epsilon (简化版本)
ax1 = axes[0, 0]
ax1.semilogx(eps_simple, abs_simple, 'bo-', linewidth=2, markersize=6)
ax1.set_xlabel(r'$\epsilon$', fontsize=12)
ax1.set_ylabel(r'$|A/B|$', fontsize=12)
ax1.set_title('Simplified: $|A/B|$ vs $\\epsilon$', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.axhline(abs_simple[-1], color='b', linestyle='--', alpha=0.5,
            label=f'ε→0 limit ≈ {abs_simple[-1]:.6f}')
ax1.legend()

# 子图2: |A/B| vs epsilon (完整版本)
ax2 = axes[0, 1]
ax2.semilogx(eps_full, abs_full, 'ro-', linewidth=2, markersize=6)
ax2.set_xlabel(r'$\epsilon$', fontsize=12)
ax2.set_ylabel(r'$|A/B|$', fontsize=12)
ax2.set_title('Full: $|A/B|$ vs $\\epsilon$', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.axhline(abs_full[-1], color='r', linestyle='--', alpha=0.5,
            label=f'ε→0 limit ≈ {abs_full[-1]:.6f}')
ax2.legend()

# 子图3: 相对于极限值的偏差 (简化版本)
ax3 = axes[1, 0]
limit_simple = abs_simple[-1]
rel_dev_simple = [(val - limit_simple) / limit_simple * 100 for val in abs_simple]
ax3.loglog(eps_simple, np.abs(rel_dev_simple), 'bs-', linewidth=2, markersize=6)
ax3.set_xlabel(r'$\epsilon$', fontsize=12)
ax3.set_ylabel('|Relative Deviation| (%)', fontsize=12)
ax3.set_title('Simplified: Convergence to ε→0 Limit', fontsize=13, fontweight='bold')
ax3.grid(True, alpha=0.3, which='both')

# 子图4: 相对于极限值的偏差 (完整版本)
ax4 = axes[1, 1]
limit_full = abs_full[-1]
rel_dev_full = [(val - limit_full) / limit_full * 100 for val in abs_full]
ax4.loglog(eps_full, np.abs(rel_dev_full), 'rs-', linewidth=2, markersize=6)
ax4.set_xlabel(r'$\epsilon$', fontsize=12)
ax4.set_ylabel('|Relative Deviation| (%)', fontsize=12)
ax4.set_title('Full: Convergence to ε→0 Limit', fontsize=13, fontweight='bold')
ax4.grid(True, alpha=0.3, which='both')

plt.suptitle(f'Epsilon Extrapolation Analysis ($\\omega={omega}$, $l={l}$)',
             fontsize=15, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])

# 保存图表
output_file = 'figures/epsilon_extrapolation.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"图表已保存: {output_file}")

# 打印详细分析
print()
print("=" * 70)
print("外推分析")
print("=" * 70)

print("\n简化版本:")
print(f"  ε = 1e-1:  |A/B| = {abs_simple[0]:.6f}")
print(f"  ε = 1e-8:  |A/B| = {abs_simple[-5]:.6f}")
print(f"  ε = 1e-9:  |A/B| = {abs_simple[-1]:.6f}")
print(f"  总变化 (1e-1 → 1e-9): {abs_simple[-1] - abs_simple[0]:.6e} ({(abs_simple[-1] - abs_simple[0])/abs_simple[0]*100:.4f}%)")
print(f"  变化 (1e-8 → 1e-9): {abs_simple[-1] - abs_simple[-5]:.6e} ({(abs_simple[-1] - abs_simple[-5])/abs_simple[-5]*100:.6f}%)")

print("\n完整版本:")
print(f"  ε = 1e-1:  |A/B| = {abs_full[0]:.6f}")
print(f"  ε = 1e-8:  |A/B| = {abs_full[-5]:.6f}")
print(f"  ε = 1e-9:  |A/B| = {abs_full[-1]:.6f}")
print(f"  总变化 (1e-1 → 1e-9): {abs_full[-1] - abs_full[0]:.6e} ({(abs_full[-1] - abs_full[0])/abs_full[0]*100:.4f}%)")
print(f"  变化 (1e-8 → 1e-9): {abs_full[-1] - abs_full[-5]:.6e} ({(abs_full[-1] - abs_full[-5])/abs_full[-5]*100:.6f}%)")

# 检查收敛性
print("\n收敛性分析:")
if len(abs_simple) >= 3:
    last_3_simple = abs_simple[-3:]
    std_simple = np.std(last_3_simple)
    print(f"  简化版本（最后3个值）: 标准差 = {std_simple:.2e} ({std_simple/np.mean(last_3_simple)*100:.6f}%)")

if len(abs_full) >= 3:
    last_3_full = abs_full[-3:]
    std_full = np.std(last_3_full)
    print(f"  完整版本（最后3个值）: 标准差 = {std_full:.2e} ({std_full/np.mean(last_3_full)*100:.6f}%)")

print("\n结论:")
if std_simple < 1e-8 and std_full < 1e-8:
    print("  ✅ 两个版本都已收敛到 ε→0 极限（标准差 < 1e-8）")
else:
    print("  ⚠️  可能需要更小的 epsilon 值来确认收敛")
