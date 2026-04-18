/**
 * @file radial_flow.hpp
 * @brief 径向流计算：Schwarzschild 背景下的标量场 (s=0)
 *
 * 通过从视界向外积分 Riccati 方程来计算反射系数 A/B。
 * 参考：径向流(check) (1).tex
 */

#ifndef RADIAL_FLOW_HPP
#define RADIAL_FLOW_HPP

#include <complex>
#include <vector>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace radial_flow {

using Complex = std::complex<double>;

/**
 * @brief Schwarzschild 黑洞参数
 */
struct SchwarzschildParams {
    double r0 = 2.0;  // Schwarzschild 半径 (2M, M=1)

    /// 度规函数 f = 1 - r0/r
    double f(double r) const {
        return 1.0 - r0 / r;
    }

    /// f 的导数：f' = r0/r^2
    double f_prime(double r) const {
        return r0 / (r * r);
    }
};

/**
 * @brief 标量波参数
 */
struct WaveParams {
    Complex omega;  // 频率（可以是复数，用于准正规模）
    int l;          // 角动量量子数
    int m = 0;      // 方位角量子数（Schwarzschild 不依赖 m）
    bool use_subleading_term = true;  // true=完整版本V=l(l+1)/r²+f'/r，false=简化版本V=l(l+1)/r²
    bool integrate_in_rstar = false;  // true=在r*坐标积分，false=在r坐标积分
};

/**
 * @brief 积分结果
 */
struct SolveResult {
    std::vector<double> r;
    std::vector<Complex> sigma;
    Complex reflection_coeff;
    bool success;
};

/**
 * @brief 径向流方程（Riccati 方程）求解器
 *
 * 变量 sigma = (J + i*omega*R) / (J - i*omega*R) 满足 Riccati ODE。
 * 在无穷远处，sigma -> -A/B 给出反射系数。
 */
class RadialFlowSolver {
public:
    RadialFlowSolver(const SchwarzschildParams& bh, const WaveParams& wave)
        : bh_(bh), wave_(wave) {}

    /**
     * @brief 势函数 V = l(l+1)/r^2 + f'/r
     *
     * 完整版本包含 f' 项，对应正确的渐近行为：
     * R ~ (A_0 + A_1/r) * r^{-1} * e^{iω r_*} + (B_0 + B_1/r) * r^{-1} * e^{-iω r_*}
     * 其中 r_* = r + r_0 ln(r/r_0 - 1) 是龟坐标
     */
    double potential_V(double r) const;

    /**
     * @brief Riccati 方程右端项（r坐标）
     * dσ/dr = -1/r - V/(2iω) + [2iω/f + V/(iω)]σ - [-1/r + V/(2iω)]σ²
     */
    Complex sigma_rhs(double r, Complex sigma) const;

    /**
     * @brief Riccati 方程右端项（r*坐标）
     * dσ/dr* = (dσ/dr) × (dr/dr*) = (dσ/dr) × (r - r₀)/r
     */
    Complex sigma_rhs_rstar(double r, Complex sigma) const;

    /**
     * @brief 计算视界附近的初始条件
     * @param epsilon 视界附近的起始偏移（相对于 r0）
     * @return (r_start, sigma_start)
     */
    std::pair<double, Complex> sigma_initial(double epsilon = 1e-6) const;

    /**
     * @brief 从视界积分到 r_max（固定步长RK4，带收敛判据）
     * @param r_max 最小积分终点（单位：r0）
     * @param epsilon 视界附近的起始偏移
     * @param n_points 输出点数
     * @param convergence_threshold 收敛阈值（震荡幅度相对于|A/B|的比例，默认0.01即1%）
     * @param n_periods 检测周期数（默认2个周期）
     * @param max_r_factor r 的最大扩展倍数（相对于 r_max，默认5倍）
     */
    SolveResult solve(double r_max = 200.0, double epsilon = 1e-6,
                      int n_points = 10000,
                      double convergence_threshold = 0.01,
                      int n_periods = 2,
                      double max_r_factor = 5.0) const;

    /**
     * @brief 从视界积分到 r_max（自适应步长RK45）
     * @param r_max 积分终点（单位：r0）
     * @param epsilon 视界附近的起始偏移
     * @param rtol 相对误差容限
     * @param atol 绝对误差容限
     */
    SolveResult solve_adaptive(double r_max = 200.0, double epsilon = 1e-6,
                               double rtol = 1e-8, double atol = 1e-10) const;

    /**
     * @brief 计算反射系数 A/B（使用自适应方法）
     */
    Complex compute_reflection_coefficient(double r_max = 200.0) const;

private:
    SchwarzschildParams bh_;
    WaveParams wave_;

    /**
     * @brief RK4 单步积分
     * @param r 当前r坐标
     * @param sigma 当前σ值（会被更新）
     * @param dr 步长（r坐标或r*坐标，取决于integrate_in_rstar；会被更新为实际的dr）
     */
    void rk4_step(double r, Complex& sigma, double& dr) const;

    /**
     * @brief RK45 单步积分（带误差估计）
     * @return (sigma_new, error_estimate)
     */
    std::pair<Complex, double> rk45_step(double r, Complex sigma, double dr) const;
};

}  // namespace radial_flow

#endif  // RADIAL_FLOW_HPP
