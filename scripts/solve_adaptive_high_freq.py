#!/usr/bin/env python3
"""
高频自适应求解器
根据频率自动选择合适的参数
"""

import numpy as np
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

def solve_adaptive_high_freq(bh, wave, r_max=200.0,
                              convergence_threshold=0.1,
                              n_periods=10,
                              max_r_factor=5.0):
    """
    自适应高频求解器

    根据频率自动选择合适的 epsilon 和 n_points

    Parameters:
    -----------
    bh : SchwarzschildParams
        黑洞参数
    wave : WaveParams
        波参数
    r_max : float
        目标最大半径（单位：r0）
    convergence_threshold : float
        收敛阈值
    n_periods : int
        检测周期数
    max_r_factor : float
        最大扩展因子

    Returns:
    --------
    result : dict
        求解结果
    """
    omega = wave.omega
    omega_abs = abs(omega)  # 使用绝对值
    solver = RadialFlowSolver(bh, wave)

    # 根据频率选择参数
    if omega_abs < 2.0:
        # 低中频：使用默认参数
        epsilon = 1e-6
        n_points = 10000
        print(f"ω = {omega_abs:.2f}: 使用标准参数 (epsilon={epsilon:.0e}, n_points={n_points})")
    elif omega_abs < 3.0:
        # 中高频：增加积分点数
        epsilon = 1e-6
        n_points = 100000
        print(f"ω = {omega_abs:.2f}: 使用增强参数 (epsilon={epsilon:.0e}, n_points={n_points})")
    elif omega_abs < 5.0:
        # 高频：增大 epsilon 和 n_points
        epsilon = 1e-3
        n_points = 400000
        print(f"ω = {omega_abs:.2f}: 使用高频参数 (epsilon={epsilon:.0e}, n_points={n_points})")
    else:
        # 极高频：使用最大参数
        epsilon = 1e-3
        n_points = 800000
        print(f"ω = {omega_abs:.2f}: 使用极高频参数 (epsilon={epsilon:.0e}, n_points={n_points})")

    # 求解
    result = solver.solve(
        r_max=r_max,
        epsilon=epsilon,
        n_points=n_points,
        convergence_threshold=convergence_threshold,
        n_periods=n_periods,
        max_r_factor=max_r_factor
    )

    return result

# 测试
if __name__ == "__main__":
    bh = SchwarzschildParams(r0=2.0)

    print("=" * 70)
    print("高频自适应求解器测试")
    print("=" * 70)
    print()

    # 测试不同频率
    test_omegas = [0.5, 1.0, 2.0, 3.0, 5.0]

    results = []

    for omega in test_omegas:
        wave = WaveParams(omega=omega, l=2)

        result = solve_adaptive_high_freq(
            bh, wave,
            r_max=200.0,
            convergence_threshold=0.1,
            n_periods=10,
            max_r_factor=5.0
        )

        r_final = result['r'][-1] / bh.r0
        converged = result['success']
        abs_reflection = np.abs(result['reflection_coeff'])

        # 检查是否有 nan
        has_nan = np.any(np.isnan(result['sigma']))

        results.append({
            'omega': omega,
            'r_final': r_final,
            'converged': converged,
            'abs_reflection': abs_reflection,
            'has_nan': has_nan
        })

        print(f"  r_final/r0 = {r_final:.1f}")
        print(f"  收敛状态 = {converged}")
        print(f"  |A/B| = {abs_reflection:.6e}")
        print(f"  包含 nan = {has_nan}")
        print()

    print("=" * 70)
    print("测试完成")
    print("=" * 70)
    print()

    # 总结
    print("总结：")
    for r in results:
        status = "✓" if (not r['has_nan'] and r['abs_reflection'] <= 1.0) else "✗"
        print(f"  ω = {r['omega']:4.1f}: {status} |A/B| = {r['abs_reflection']:.6e}, "
              f"converged = {r['converged']}, nan = {r['has_nan']}")
