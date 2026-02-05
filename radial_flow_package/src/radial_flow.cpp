/**
 * @file radial_flow.cpp
 * @brief 径向流计算实现
 */

#include "radial_flow.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace radial_flow {

double RadialFlowSolver::potential_V(double r) const {
    int l = wave_.l;
    double V_angular = static_cast<double>(l * (l + 1)) / (r * r);

    if (wave_.use_subleading_term) {
        // 完整版本：V = l(l+1)/r^2 + f'/r
        // 对应正确的渐近行为：R ~ (A_0 + A_1/r) * r^{-1} * e^{iω r_*}+ (B_0 + B_1/r) * r^{-1} * e^{-iω r_*}
        double f_prime = bh_.f_prime(r);
        return V_angular + f_prime / r;
    } else {
        // 简化版本：V = l(l+1)/r^2
        // 对应简化的渐近行为：R ~ A * e^{iω r} + B * e^{-iω r}
        return V_angular;
    }
}

Complex RadialFlowSolver::sigma_rhs(double r, Complex sigma) const {
    Complex omega = wave_.omega;
    double f = bh_.f(r);
    double V = potential_V(r);

    // 系数（根据势函数版本决定是否包含 -1/r 项）
    Complex a, b, c;

    if (wave_.use_subleading_term) {
        // 完整版本：去掉 a 系数的 -1/r 项，c 系数保持负号
        a = - V/(2.0*Complex(0, 1)*omega);
        b = 2.0*Complex(0, 1)*omega/f + V/(Complex(0, 1)*omega);
        c = -(-1.0/r + V/(2.0*Complex(0, 1)*omega));
    } else {
        // 简化版本：保留 -1/r 项
        a = -1.0/r - V/(2.0*Complex(0, 1)*omega);
        b = 2.0*Complex(0, 1)*omega/f + V/(Complex(0, 1)*omega);
        c = -(-1.0/r + V/(2.0*Complex(0, 1)*omega));
    }

    return a + b*sigma + c*sigma*sigma;
}

std::pair<double, Complex> RadialFlowSolver::sigma_initial(double epsilon) const {
    double r0 = bh_.r0;
    Complex omega = wave_.omega;
    int l = wave_.l;

    // 视界处：σ0 = 0（手动设置，确保纯入射波）
    Complex sigma0 = 0.0;

    // σ0' 的表达式（根据势函数版本选择不同的边界条件）
    double f0_prime = 1.0 / (r0);
    Complex sigma0_prime;

    if (wave_.use_subleading_term) {
        // 完整版本（去掉 -1/r 项后的新边界条件）
        // σ₀' = [r₀⁻²·l(l+1) + r₀⁻¹·f₀'] · f₀' / [(2iω) · (2iω - f₀')]
        //     = V₀ · f₀' / [(2iω) · (2iω - f₀')]
        // 其中 V₀ = l(l+1)/r₀² + f₀'/r₀ 是视界处的势函数值
        double V0 = static_cast<double>(l*(l+1))/(r0*r0) + f0_prime/r0;
        Complex numerator = V0 * f0_prime;
        Complex denominator = 2.0*Complex(0,1)*omega * (2.0*Complex(0,1)*omega - f0_prime);
        sigma0_prime = numerator / denominator;
    } else {
        // 简化版本（保留原来的边界条件，.tex 第105行）
        Complex numerator = -f0_prime * (1.0 + static_cast<double>(l*(l+1))/(2.0*Complex(0,1)*omega*r0));
        Complex denominator = r0*f0_prime - 2.0*Complex(0,1)*omega*r0;
        sigma0_prime = numerator / denominator;
    }

    double r_start = r0 + epsilon * r0;
    Complex sigma_start = sigma0 + sigma0_prime * epsilon * r0;

    return {r_start, sigma_start};
}

void RadialFlowSolver::rk4_step(double r, Complex& sigma, double dr) const {
    Complex k1 = sigma_rhs(r, sigma);
    Complex k2 = sigma_rhs(r + 0.5*dr, sigma + 0.5*dr*k1);
    Complex k3 = sigma_rhs(r + 0.5*dr, sigma + 0.5*dr*k2);
    Complex k4 = sigma_rhs(r + dr, sigma + dr*k3);

    sigma = sigma + (dr/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

SolveResult RadialFlowSolver::solve(double r_max, double epsilon, int n_points,
                                     double convergence_threshold, int n_periods,
                                     double max_r_factor) const {
    auto [r_start, sigma_start] = sigma_initial(epsilon);

    double r_end = r_max * bh_.r0;
    double dr = (r_end - r_start) / static_cast<double>(n_points - 1);

    std::vector<double> r_vec;
    std::vector<Complex> sigma_vec;

    r_vec.reserve(n_points * 2);
    sigma_vec.reserve(n_points * 2);

    double r = r_start;
    Complex sigma = sigma_start;

    // 计算检测窗口大小
    double omega_abs = std::abs(wave_.omega);
    double period = 2.0 * M_PI / omega_abs;
    int points_per_period = static_cast<int>(period / dr) + 1;
    int check_window = n_periods * points_per_period;

    // 调试输出（已禁用）
    // std::cout << "Debug: omega=" << omega_abs << ", period=" << period
    //           << ", points_per_period=" << points_per_period
    //           << ", check_window=" << check_window << std::endl;

    double max_r_end = max_r_factor * r_end;
    bool converged = false;
    int check_count = 0;  // 检测次数计数器

    // 积分并检测收敛
    while (r < max_r_end && !converged) {
        r_vec.push_back(r);
        sigma_vec.push_back(sigma);

        // 当有足够的点时开始检测收敛（不限制必须 r > r_max）
        if (static_cast<int>(sigma_vec.size()) > check_window) {
            check_count++;

            // 提取最近 check_window 个点
            int start_idx = sigma_vec.size() - check_window;

            // 应用相位修正后再检测收敛
            std::vector<Complex> recent_sigma_corrected;
            recent_sigma_corrected.reserve(check_window);

            for (int i = start_idx; i < static_cast<int>(sigma_vec.size()); ++i) {
                double r_i = r_vec[i];
                double r_star = r_i + bh_.r0 * std::log(r_i / bh_.r0 - 1.0);
                Complex phase_correction = std::exp(Complex(0, -2.0) * wave_.omega * r_star);
                recent_sigma_corrected.push_back(sigma_vec[i] * phase_correction);
            }

            // 计算平均值（使用修正后的 sigma）
            Complex sigma_mean = 0.0;
            for (const auto& s : recent_sigma_corrected) {
                sigma_mean += s;
            }
            sigma_mean /= static_cast<double>(recent_sigma_corrected.size());

            // 计算震荡幅度
            double max_deviation = 0.0;
            for (const auto& s : recent_sigma_corrected) {
                double deviation = std::abs(s - sigma_mean);
                if (deviation > max_deviation) {
                    max_deviation = deviation;
                }
            }

            // 计算相对震荡幅度
            double sigma_mean_abs = std::abs(sigma_mean);
            if (sigma_mean_abs > 1e-10) {
                double relative_oscillation = max_deviation / sigma_mean_abs;

                // 每1000次检测输出一次（已禁用）
                // if (check_count % 1000 == 0) {
                //     std::cout << "Debug: check_count=" << check_count
                //               << ", r/r0=" << r/bh_.r0
                //               << ", relative_osc=" << relative_oscillation
                //               << ", threshold=" << convergence_threshold
                //               << ", r>r_end=" << (r > r_end) << std::endl;
                // }

                // 判断是否收敛（但至少要积分到 r_max）
                if (r > r_end && relative_oscillation < convergence_threshold) {
                    // std::cout << "Converged! r/r0=" << r/bh_.r0
                    //           << ", relative_osc=" << relative_oscillation << std::endl;
                    converged = true;
                    break;  // 收敛后立即停止
                }
            }
        }

        // 继续积分
        rk4_step(r, sigma, dr);
        r += dr;
    }

    // 添加最后一个点
    if (!converged || r_vec.empty() || r_vec.back() != r) {
        r_vec.push_back(r);
        sigma_vec.push_back(sigma);
    }

    // 反射系数 = -sigma (在无穷远处)
    Complex reflection_coeff = -sigma_vec.back();

    return {r_vec, sigma_vec, reflection_coeff, converged};
}

Complex RadialFlowSolver::compute_reflection_coefficient(double r_max) const {
    SolveResult result = solve(r_max);
    return result.reflection_coeff;
}

}  // namespace radial_flow
