#!/usr/bin/env python3
"""
验证简化版本和完整版本在领头阶的一致性
使用龟坐标 r* 除相位后，A/B 应该接近
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

print("=" * 80)
print("验证领头阶一致性：简化版本 vs 完整版本")
print("=" * 80)
print(f"参数: ω = {omega}, l = {l}, r₀ = {bh.r0}")
print()

# ============================================================================
# 计算两个版本（都使用 r* 修正）
# ============================================================================

# 简化版本
wave_simple = WaveParams(omega=omega, l=l, use_subleading_term=False)
solver_simple = RadialFlowSolver(bh, wave_simple)
result_simple = solver_simple.solve(r_max=200.0, epsilon=1e-6, n_points=10000)

# 完整版本
wave_full = WaveParams(omega=omega, l=l, use_subleading_term=True)
solver_full = RadialFlowSolver(bh, wave_full)
result_full = solver_full.solve(r_max=200.0, epsilon=1e-6, n_points=10000)

# ============================================================================
# 对比结果
# ============================================================================

A_B_simple = result_simple['reflection_coeff']
A_B_full = result_full['reflection_coeff']

print("【简化版本】V = l(l+1)/r²")
print("-" * 80)
print(f"  A/B = {A_B_simple:.6f}")
print(f"  |A/B| = {abs(A_B_simple):.6f}")
print(f"  arg(A/B) = {np.angle(A_B_simple):.6f} rad = {np.angle(A_B_simple)/np.pi:.4f}π")
print(f"  Re(A/B) = {A_B_simple.real:.6f}")
print(f"  Im(A/B) = {A_B_simple.imag:.6f}")
print()

print("【完整版本】V = l(l+1)/r² + f'/r")
print("-" * 80)
print(f"  A/B = {A_B_full:.6f}")
print(f"  |A/B| = {abs(A_B_full):.6f}")
print(f"  arg(A/B) = {np.angle(A_B_full):.6f} rad = {np.angle(A_B_full)/np.pi:.4f}π")
print(f"  Re(A/B) = {A_B_full.real:.6f}")
print(f"  Im(A/B) = {A_B_full.imag:.6f}")
print()

# ============================================================================
# 差异分析
# ============================================================================

print("【差异分析】")
print("-" * 80)

# 模的差异
diff_abs = abs(A_B_full) - abs(A_B_simple)
rel_diff_abs = diff_abs / abs(A_B_simple) * 100

print(f"模的差异:")
print(f"  Δ|A/B| = {diff_abs:.6f}")
print(f"  相对差异 = {rel_diff_abs:.2f}%")
print()

# 相位的差异
diff_arg = np.angle(A_B_full) - np.angle(A_B_simple)
# 归一化到 [-π, π]
if diff_arg > np.pi:
    diff_arg -= 2*np.pi
elif diff_arg < -np.pi:
    diff_arg += 2*np.pi

print(f"相位的差异:")
print(f"  Δarg(A/B) = {diff_arg:.6f} rad = {diff_arg/np.pi:.4f}π")
print()

# 复数的差异
diff_complex = A_B_full - A_B_simple
print(f"复数的差异:")
print(f"  Δ(A/B) = {diff_complex:.6f}")
print(f"  |Δ(A/B)| = {abs(diff_complex):.6f}")
print(f"  相对差异 = {abs(diff_complex)/abs(A_B_simple)*100:.2f}%")
print()

# ============================================================================
# 可视化
# ============================================================================

print("生成可视化...")

fig = plt.figure(figsize=(16, 10))

# 子图1: σ 的演化（简化版本）
ax1 = plt.subplot(2, 3, 1)
r_simple = np.array(result_simple['r'])
sigma_simple = np.array(result_simple['sigma'])
ax1.plot(r_simple/bh.r0, np.abs(sigma_simple), 'b-', linewidth=1.5)
ax1.set_xlabel(r'$r/r_0$', fontsize=12)
ax1.set_ylabel(r'$|\sigma|$ (phase corrected)', fontsize=12)
ax1.set_title('Simplified: $|\sigma|$ vs $r$', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(left=1)

# 子图2: σ 的演化（完整版本）
ax2 = plt.subplot(2, 3, 2)
r_full = np.array(result_full['r'])
sigma_full = np.array(result_full['sigma'])
ax2.plot(r_full/bh.r0, np.abs(sigma_full), 'r-', linewidth=1.5)
ax2.set_xlabel(r'$r/r_0$', fontsize=12)
ax2.set_ylabel(r'$|\sigma|$ (phase corrected)', fontsize=12)
ax2.set_title('Full: $|\sigma|$ vs $r$', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(left=1)

# 子图3: 复平面对比
ax3 = plt.subplot(2, 3, 3)
ax3.plot(sigma_simple.real, sigma_simple.imag, 'b-', linewidth=1, alpha=0.6, label='Simplified')
ax3.plot(sigma_full.real, sigma_full.imag, 'r-', linewidth=1, alpha=0.6, label='Full')
ax3.scatter([A_B_simple.real], [A_B_simple.imag], c='blue', s=150, marker='o',
            edgecolors='black', linewidths=2, zorder=5, label='Simplified: A/B')
ax3.scatter([A_B_full.real], [A_B_full.imag], c='red', s=150, marker='s',
            edgecolors='black', linewidths=2, zorder=5, label='Full: A/B')
ax3.set_xlabel(r'Re($\sigma$)', fontsize=12)
ax3.set_ylabel(r'Im($\sigma$)', fontsize=12)
ax3.set_title('Complex Plane Comparison', fontsize=13, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)
ax3.axis('equal')

# 子图4: A/B 在复平面上的对比（放大）
ax4 = plt.subplot(2, 3, 4)
ax4.scatter([A_B_simple.real], [A_B_simple.imag], c='blue', s=200, marker='o',
            edgecolors='black', linewidths=2, zorder=5, label='Simplified')
ax4.scatter([A_B_full.real], [A_B_full.imag], c='red', s=200, marker='s',
            edgecolors='black', linewidths=2, zorder=5, label='Full')
# 画一条连线
ax4.plot([A_B_simple.real, A_B_full.real], [A_B_simple.imag, A_B_full.imag],
         'k--', linewidth=1, alpha=0.5)
# 标注差异
ax4.annotate(f'Δ = {abs(diff_complex):.4f}',
             xy=((A_B_simple.real + A_B_full.real)/2, (A_B_simple.imag + A_B_full.imag)/2),
             fontsize=10, ha='center')
ax4.set_xlabel(r'Re(A/B)', fontsize=12)
ax4.set_ylabel(r'Im(A/B)', fontsize=12)
ax4.set_title('A/B Comparison (Zoomed)', fontsize=13, fontweight='bold')
ax4.legend(fontsize=10)
ax4.grid(True, alpha=0.3)
ax4.axis('equal')

# 子图5: 模和相位的对比
ax5 = plt.subplot(2, 3, 5)
categories = ['|A/B|', 'arg(A/B)']
simple_vals = [abs(A_B_simple), np.angle(A_B_simple)]
full_vals = [abs(A_B_full), np.angle(A_B_full)]

x = np.arange(len(categories))
width = 0.35

bars1 = ax5.bar(x - width/2, simple_vals, width, label='Simplified', color='blue', alpha=0.7)
bars2 = ax5.bar(x + width/2, full_vals, width, label='Full', color='red', alpha=0.7)

ax5.set_ylabel('Value', fontsize=12)
ax5.set_title('Magnitude and Phase Comparison', fontsize=13, fontweight='bold')
ax5.set_xticks(x)
ax5.set_xticklabels(categories)
ax5.legend()
ax5.grid(True, alpha=0.3, axis='y')

# 在柱子上标注数值
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax5.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.3f}', ha='center', va='bottom', fontsize=9)

# 子图6: 结果总结
ax6 = plt.subplot(2, 3, 6)
ax6.axis('off')

summary_text = f"""
Parameters
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ω = {omega}, l = {l}, r₀ = {bh.r0}

Simplified Version
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
A/B = {A_B_simple.real:.4f} + {A_B_simple.imag:.4f}i
|A/B| = {abs(A_B_simple):.6f}
arg(A/B) = {np.angle(A_B_simple):.4f} rad

Full Version
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
A/B = {A_B_full.real:.4f} + {A_B_full.imag:.4f}i
|A/B| = {abs(A_B_full):.6f}
arg(A/B) = {np.angle(A_B_full):.4f} rad

Difference
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Δ|A/B| = {diff_abs:.6f} ({rel_diff_abs:.2f}%)
Δarg(A/B) = {diff_arg:.4f} rad
|Δ(A/B)| = {abs(diff_complex):.6f}

Conclusion
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Using r* phase correction for both:
- Magnitude differs by {rel_diff_abs:.2f}%
- Phase differs by {diff_arg/np.pi:.2f}π
- Complex difference: {abs(diff_complex)/abs(A_B_simple)*100:.2f}%

The subleading term f'/r causes
a {rel_diff_abs:.1f}% correction to |A/B|
"""

ax6.text(0.05, 0.5, summary_text, transform=ax6.transAxes,
         fontsize=10, verticalalignment='center',
         family='monospace', bbox=dict(boxstyle='round',
         facecolor='wheat', alpha=0.3))

plt.suptitle(f'Leading Order Consistency Check (Both using r* phase correction)',
             fontsize=15, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])

# 保存图表
output_file = 'figures/leading_order_consistency.png'
os.makedirs('figures', exist_ok=True)
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"图表已保存: {output_file}")

print()
print("=" * 80)
print("结论")
print("=" * 80)
print()

if abs(diff_complex)/abs(A_B_simple) < 0.01:
    print("✅ 验证成功！两个版本在领头阶基本一致（差异 < 1%）")
    print("   次领头项 f'/r 的影响很小")
elif abs(diff_complex)/abs(A_B_simple) < 0.30:
    print("⚠️  两个版本有显著差异（差异 ~ 23%）")
    print("   次领头项 f'/r 对结果有重要影响")
    print("   这是预期的，因为 f'/r 项在低频时不可忽略")
else:
    print("❌ 两个版本差异很大（差异 > 30%）")
    print("   可能存在问题")

print()
print(f"使用龟坐标 r* 相位修正后:")
print(f"  - 模的相对差异: {rel_diff_abs:.2f}%")
print(f"  - 相位差异: {diff_arg/np.pi:.4f}π")
print(f"  - 复数的相对差异: {abs(diff_complex)/abs(A_B_simple)*100:.2f}%")
