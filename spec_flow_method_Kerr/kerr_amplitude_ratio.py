#!/usr/bin/env python3
"""
Kerr黑洞引力波散射振幅比计算器
输入: a (自旋), omega (频率), l (角量子数), m (方位角量子数)
输出: B_inc/B_ref (复数), lambda (分离常数)
"""
import sys
sys.path.insert(0, 'kerr_matcher_project/src')

from kerr_matcher.params import SolverParams
from kerr_matcher.solver import solve_case


def compute_amplitude_ratio(a, omega, l, m, lambda_sep=None, r_match=8.0, n_cheb=32):
    """
    计算Kerr黑洞引力波散射振幅比

    参数:
        a: 自旋参数 (0 <= a < M, M=1)
        omega: 频率
        l: 角量子数 (l >= 2)
        m: 方位角量子数 (|m| <= l)
        lambda_sep: 球旋分离常数 (None则自动计算)
        r_match: 匹配半径
        n_cheb: Chebyshev多项式阶数

    返回:
        dict: {
            'ratio': complex,  # B_inc/B_ref (入射/反射)
            'lambda': complex,  # 分离常数
            'ratio_abs': float,  # |ratio|
            'ratio_arg': float   # arg(ratio)
        }
    """
    params = SolverParams(
        M=1.0, a=a, omega=omega, ell=l, m_mode=m,
        lambda_sep=lambda_sep, r_match=r_match, n_cheb=n_cheb, flow_eps=1e-6
    )

    result = solve_case(params)

    # ratio_up_over_um = R_+/R_- = 反射/入射
    # 所以 B_inc/B_ref = 1/ratio_up_over_um
    ratio_ref_over_inc = result.spectral.ratio_up_over_um
    ratio_inc_over_ref = 1.0 / ratio_ref_over_inc

    return {
        'ratio': ratio_inc_over_ref,
        'lambda': params.lambda_value,
        'ratio_abs': abs(ratio_inc_over_ref),
        'ratio_arg': ratio_inc_over_ref.real / abs(ratio_inc_over_ref) if abs(ratio_inc_over_ref) > 0 else 0
    }


if __name__ == '__main__':
    # 示例调用
    if len(sys.argv) >= 5:
        a = float(sys.argv[1])
        omega = float(sys.argv[2])
        l = int(sys.argv[3])
        m = int(sys.argv[4])
        lambda_sep = float(sys.argv[5]) if len(sys.argv) > 5 else None
    else:
        # 默认参数
        a, omega, l, m = 0.3, 0.2, 2, 2
        lambda_sep = 3.601

    result = compute_amplitude_ratio(a, omega, l, m, lambda_sep)

    print(f"输入参数: a={a}, ω={omega}, l={l}, m={m}")
    print(f"λ = {result['lambda']}")
    print(f"B_inc/B_ref = {result['ratio']}")
    print(f"|B_inc/B_ref| = {result['ratio_abs']:.6e}")
