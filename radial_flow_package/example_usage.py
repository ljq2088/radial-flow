#!/usr/bin/env python3
"""
径向流计算 - 简单使用示例
"""

import sys
import os

# 添加当前目录到 Python 路径（使用相对路径）
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver
import numpy as np

# ============================================================================
# 基本用法
# ============================================================================

# 1. 设置黑洞参数
bh = SchwarzschildParams(r0=2.0)  # Schwarzschild 半径 (M=1 时 r0=2M=2)

# 2. 设置波参数
omega = 0.2  # 频率
l = 0        # 角动量量子数

# 3. 选择势函数版本
# 简化版本: V = l(l+1)/r² (保留 -1/r 项)
wave_simple = WaveParams(omega=omega, l=l, use_subleading_term=False)

# 完整版本: V = l(l+1)/r² + f'/r (去掉 a 系数的 -1/r 项)
wave_full = WaveParams(omega=omega, l=l, use_subleading_term=True)

# 4. 创建求解器并计算
solver = RadialFlowSolver(bh, wave_full)  # 使用完整版本

result = solver.solve(
    r_max=200.0,                    # 最小积分终点 (单位: r0)
    epsilon=1e-6,                   # 视界偏移参数
    n_points=10000,                 # 输出点数
    convergence_threshold=0.01,     # 收敛阈值 (1%)
    n_periods=2,                    # 检测周期数
    max_r_factor=2.0                # 最大扩展倍数
)

# 5. 查看结果
print(f"收敛状态: {result['success']}")
print(f"反射系数 A/B = {result['reflection_coeff']:.6f}")
print(f"|A/B| = {abs(result['reflection_coeff']):.6f}")
print(f"相位 arg(A/B) = {np.angle(result['reflection_coeff']):.6f} rad")
print(f"实际 r_max = {result['r'][-1]:.2f}")

# ============================================================================
# 高级用法：频率扫描
# ============================================================================

def frequency_scan(l=2, omegas=None):
    """扫描不同频率的反射系数"""
    if omegas is None:
        omegas = np.linspace(0.1, 2.0, 20)

    bh = SchwarzschildParams(r0=2.0)
    results = []

    print(f"\n频率扫描 (l={l}):")
    print("-" * 60)

    for omega in omegas:
        wave = WaveParams(omega=omega, l=l, use_subleading_term=True)
        solver = RadialFlowSolver(bh, wave)
        result = solver.solve(r_max=200.0)

        A_over_B = result['reflection_coeff']
        results.append({
            'omega': omega,
            'A_over_B': A_over_B,
            'abs': abs(A_over_B),
            'arg': np.angle(A_over_B)
        })

        print(f"ω={omega:.2f}: |A/B|={abs(A_over_B):.4f}, arg={np.angle(A_over_B):.4f}")

    return results

# ============================================================================
# 参数说明
# ============================================================================
"""
WaveParams 参数:
  - omega: 频率 (复数或实数)
  - l: 角动量量子数 (整数, l >= 0)
  - m: 方位角量子数 (默认 0)
  - use_subleading_term:
      * True  = 完整版本 V = l(l+1)/r² + f'/r (去掉 a 系数的 -1/r)
      * False = 简化版本 V = l(l+1)/r² (保留 -1/r 项)

solve() 参数:
  - r_max: 最小积分终点 (单位: r0)
  - epsilon: 视界偏移 (默认 1e-6, 高频时可用 1e-3)
  - n_points: 输出点数
  - convergence_threshold: 收敛阈值 (相对震荡幅度)
  - n_periods: 检测周期数
  - max_r_factor: 最大扩展倍数

返回结果:
  - r: 半径数组
  - sigma: 修正后的 sigma 数组 (已应用相位修正)
  - reflection_coeff: 反射系数 A/B
  - success: 是否成功收敛
"""

if __name__ == "__main__":
    # 运行示例
    print("径向流计算示例")
    print("=" * 80)

    # 单点计算
    bh = SchwarzschildParams(r0=2.0)
    wave = WaveParams(omega=0.2, l=0, use_subleading_term=True)
    solver = RadialFlowSolver(bh, wave)
    result = solver.solve(r_max=200.0)

    print(f"\n单点计算 (ω={wave.omega}, l={wave.l}):")
    print(f"  |A/B| = {abs(result['reflection_coeff']):.6f}")

    # 频率扫描示例（取消注释以运行）
    # frequency_scan(l=2)
