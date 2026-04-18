#!/usr/bin/env python3
"""
分析为什么无穷远处定义的 A/B 会在方程上出现差异
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
print("分析：为什么无穷远处定义的 A/B 会在方程上出现差异")
print("=" * 80)
print()

# ============================================================================
# 1. 势函数的差异
# ============================================================================

print("【1. 势函数的差异】")
print("-" * 80)

r_array = np.logspace(np.log10(2.1), np.log10(1000), 1000)

V_simple = l*(l+1) / r_array**2
V_full = l*(l+1) / r_array**2 + (bh.r0/r_array**2) / r_array

print(f"在不同距离处的势函数值：")
print()
for r_test in [5, 10, 50, 100, 500]:
    V_s = l*(l+1) / r_test**2
    V_f = l*(l+1) / r_test**2 + (bh.r0/r_test**2) / r_test
    diff = V_f - V_s
    rel_diff = diff / V_s * 100 if V_s > 0 else 0
    print(f"  r = {r_test:4.0f} r₀: V_simple = {V_s:.6e}, V_full = {V_f:.6e}")
    print(f"           差异 = {diff:.6e} ({rel_diff:+.1f}%)")
    print()

# ============================================================================
# 2. 积分路径的贡献
# ============================================================================

print("【2. 积分路径的贡献】")
print("-" * 80)

# 计算两个版本
wave_simple = WaveParams(omega=omega, l=l, use_subleading_term=False)
solver_simple = RadialFlowSolver(bh, wave_simple)
result_simple = solver_simple._cpp_solver.solve(200.0, 1e-6, 10000, 0.01, 2, 2.0)

wave_full = WaveParams(omega=omega, l=l, use_subleading_term=True)
solver_full = RadialFlowSolver(bh, wave_full)
result_full = solver_full._cpp_solver.solve(200.0, 1e-6, 10000, 0.01, 2, 2.0)

r_simple = np.array(result_simple.r)
sigma_simple_raw = np.array(result_simple.sigma)

r_full = np.array(result_full.r)
sigma_full_raw = np.array(result_full.sigma)

print(f"简化版本:")
print(f"  积分起点: r = {r_simple[0]:.2f}")
print(f"  积分终点: r = {r_simple[-1]:.2f}")
print(f"  积分点数: {len(r_simple)}")
print()

print(f"完整版本:")
print(f"  积分起点: r = {r_full[0]:.2f}")
print(f"  积分终点: r = {r_full[-1]:.2f}")
print(f"  积分点数: {len(r_full)}")
print()

# ============================================================================
# 3. |σ| 的演化对比
# ============================================================================

print("【3. |σ| 的演化对比】")
print("-" * 80)

# 在几个关键点对比 |σ|
test_points = [10, 50, 100, 200, 400]

print(f"在不同距离处的 |σ| 值：")
print()

for r_test in test_points:
    # 找到最接近的点
    idx_simple = np.argmin(np.abs(r_simple - r_test))
    idx_full = np.argmin(np.abs(r_full - r_test))

    if r_simple[idx_simple] < r_test * 1.1 and r_full[idx_full] < r_test * 1.1:
        sigma_s = abs(sigma_simple_raw[idx_simple])
        sigma_f = abs(sigma_full_raw[idx_full])
        diff = sigma_f - sigma_s
        rel_diff = diff / sigma_s * 100

        print(f"  r ≈ {r_test:3.0f}: |σ_simple| = {sigma_s:.6f}, |σ_full| = {sigma_f:.6f}")
        print(f"           差异 = {diff:+.6f} ({rel_diff:+.2f}%)")
        print()

# ============================================================================
# 4. 物理解释
# ============================================================================

print("=" * 80)
print("物理解释")
print("=" * 80)
print()

print("关键点：")
print()
print("1. A/B 虽然在无穷远处定义，但它的值是通过从视界到无穷远的")
print("   整个积分路径决定的")
print()
print("2. Riccati 方程：dσ/dr = a + bσ + cσ²")
print("   - 系数 a, b, c 依赖于势函数 V(r)")
print("   - 简化版本：V = l(l+1)/r²")
print("   - 完整版本：V = l(l+1)/r² + f'/r")
print()
print("3. f'/r 项在整个积分路径上都有贡献：")
print("   - 在 r ~ 10 r₀ 时，f'/r ~ r₀/r³ 仍然显著")
print("   - 这会影响 σ 的演化")
print("   - 累积效应导致最终的 |σ| 不同")
print()
print("4. 相位修正只是在最后除掉一个相位因子：")
print("   - σ_corrected = σ_raw × e^{-2iωr*}")
print("   - 这只改变相位，不改变模")
print("   - 但 |σ_raw| 已经在积分过程中产生了差异")
print()
print("5. 为什么相位接近但模不同？")
print("   - 相位主要由远场渐近行为决定")
print("   - 在 r → ∞ 时，V → 0，两个版本的方程趋于相同")
print("   - 因此相位演化相似")
print("   - 但模受整个积分路径影响")
print("   - f'/r 在中等距离的贡献累积起来，导致模的差异")
print()

# ============================================================================
# 5. 可视化
# ============================================================================

print("生成可视化...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 子图1: 势函数对比
ax1 = axes[0, 0]
ax1.loglog(r_array/bh.r0, V_simple, 'b-', linewidth=2, label='Simplified: $l(l+1)/r^2$')
ax1.loglog(r_array/bh.r0, V_full, 'r-', linewidth=2, label="Full: $l(l+1)/r^2 + f'/r$")
ax1.loglog(r_array/bh.r0, np.abs(V_full - V_simple), 'g--', linewidth=2, label='Difference')
ax1.set_xlabel(r'$r/r_0$', fontsize=12)
ax1.set_ylabel(r'$V(r)$', fontsize=12)
ax1.set_title('Potential Comparison', fontsize=13, fontweight='bold')
ax1.legend()
ax1.grid(True, alpha=0.3, which='both')

# 子图2: |σ| 的演化对比
ax2 = axes[0, 1]
ax2.plot(r_simple/bh.r0, np.abs(sigma_simple_raw), 'b-', linewidth=1.5, label='Simplified', alpha=0.7)
ax2.plot(r_full/bh.r0, np.abs(sigma_full_raw), 'r-', linewidth=1.5, label='Full', alpha=0.7)
ax2.set_xlabel(r'$r/r_0$', fontsize=12)
ax2.set_ylabel(r'$|\sigma|$ (raw, before phase correction)', fontsize=12)
ax2.set_title(r'$|\sigma|$ Evolution', fontsize=13, fontweight='bold')
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xlim(left=1)

# 子图3: |σ| 的相对差异
ax3 = axes[1, 0]
# 插值到相同的 r 网格
r_common = np.linspace(max(r_simple[0], r_full[0]), min(r_simple[-1], r_full[-1]), 1000)
sigma_simple_interp = np.interp(r_common, r_simple, np.abs(sigma_simple_raw))
sigma_full_interp = np.interp(r_common, r_full, np.abs(sigma_full_raw))
rel_diff = (sigma_full_interp - sigma_simple_interp) / sigma_simple_interp * 100

ax3.plot(r_common/bh.r0, rel_diff, 'g-', linewidth=2)
ax3.axhline(0, color='k', linestyle='--', alpha=0.3)
ax3.set_xlabel(r'$r/r_0$', fontsize=12)
ax3.set_ylabel(r'Relative Difference (%)', fontsize=12)
ax3.set_title(r'$|\sigma|$ Relative Difference', fontsize=13, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.set_xlim(left=1)

# 子图4: 说明文字
ax4 = axes[1, 1]
ax4.axis('off')

explanation = """
Why does A/B differ even though
it's defined at infinity?

Key Points:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. A/B is defined at r → ∞, but
   its value is determined by the
   ENTIRE integration path from
   the horizon to infinity

2. The subleading term f'/r has
   contributions throughout the
   integration path, especially
   at intermediate distances
   (r ~ 10-100 r₀)

3. These contributions accumulate
   and affect |σ|, which in turn
   affects |A/B|

4. Phase correction only removes
   a phase factor at the end:
   σ_corrected = σ_raw × e^(-2iωr*)

   This doesn't change |σ_raw|,
   which already differs due to
   the different integration paths

5. Phase is similar because it's
   mainly determined by far-field
   behavior where V → 0 for both

6. Magnitude differs because it's
   affected by the entire path
   where f'/r is non-negligible
"""

ax4.text(0.05, 0.5, explanation, transform=ax4.transAxes,
         fontsize=10, verticalalignment='center',
         family='monospace', bbox=dict(boxstyle='round',
         facecolor='lightblue', alpha=0.3))

plt.suptitle('Why A/B Differs: Integration Path vs Asymptotic Definition',
             fontsize=14, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])

output_file = 'figures/asymptotic_vs_integration_path.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"图表已保存: {output_file}")

print()
print("=" * 80)
print("总结")
print("=" * 80)
print()
print("A/B 虽然在无穷远处定义，但它的值是通过整个积分路径决定的。")
print("f'/r 项在整个路径上的累积贡献导致了 |A/B| 的差异。")
print("这不是矛盾，而是积分方程的本质特征。")
