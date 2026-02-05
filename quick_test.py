#!/usr/bin/env python3
"""
快速测试脚本 - 演示如何调用径向流计算
"""

import sys
sys.path.insert(0, '/home/ljq/code/radial_flow')

from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver
import numpy as np

# ============================================================================
# 基本调用示例
# ============================================================================

# 步骤1: 设置黑洞参数
bh = SchwarzschildParams(r0=2.0)  # Schwarzschild半径，M=1时r0=2M=2

# 步骤2: 设置波参数
omega = 0.2  # 频率
l = 0        # 角动量量子数

# 步骤3: 选择是否包含次领头项 (subleading term: f'/r)
# False = 简化版本: V = l(l+1)/r²
# True  = 完整版本: V = l(l+1)/r² + f'/r
wave = WaveParams(omega=omega, l=l, use_subleading_term=True)

# 步骤4: 创建求解器
solver = RadialFlowSolver(bh, wave)

# 步骤5: 求解（使用默认参数）
result = solver.solve()

# 步骤6: 查看结果
print(f"反射系数 A/B = {result['reflection_coeff']:.6f}")
print(f"|A/B| = {abs(result['reflection_coeff']):.6f}")
print(f"相位 arg(A/B) = {np.angle(result['reflection_coeff']):.6f} rad")
print(f"收敛状态: {result['success']}")
print(f"积分终点: r_max = {result['r'][-1]:.2f}")

# ============================================================================
# 自定义参数调用
# ============================================================================

print("\n" + "="*60)
print("自定义参数示例")
print("="*60)

result_custom = solver.solve(
    r_max=200.0,                    # 最小积分终点（单位: r0）
    epsilon=1e-6,                   # 视界偏移参数
    n_points=10000,                 # 输出点数
    convergence_threshold=0.01,     # 收敛阈值（1%）
    n_periods=2,                    # 检测周期数
    max_r_factor=2.0                # 最大扩展倍数
)

print(f"|A/B| = {abs(result_custom['reflection_coeff']):.6f}")

# ============================================================================
# 对比简化版本和完整版本
# ============================================================================

print("\n" + "="*60)
print("简化版本 vs 完整版本对比")
print("="*60)

# 简化版本
wave_simple = WaveParams(omega=0.2, l=0, use_subleading_term=False)
solver_simple = RadialFlowSolver(bh, wave_simple)
result_simple = solver_simple.solve(r_max=200.0, n_points=10000)

# 完整版本
wave_full = WaveParams(omega=0.2, l=0, use_subleading_term=True)
solver_full = RadialFlowSolver(bh, wave_full)
result_full = solver_full.solve(r_max=200.0, n_points=10000)

print(f"简化版本: |A/B| = {abs(result_simple['reflection_coeff']):.6f}")
print(f"完整版本: |A/B| = {abs(result_full['reflection_coeff']):.6f}")
diff = abs(abs(result_simple['reflection_coeff']) - abs(result_full['reflection_coeff']))
print(f"差异: {diff:.6f} ({diff/abs(result_simple['reflection_coeff'])*100:.2f}%)")
