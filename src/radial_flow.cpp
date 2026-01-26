/**
 * @file radial_flow.cpp
 * @brief 径向流计算实现
 */

#include "radial_flow.hpp"
#include <iostream>
#include <iomanip>

namespace radial_flow {

double RadialFlowSolver::potential_V(double r) const {
    int l = wave_.l;
    return static_cast<double>(l * (l + 1)) / (r * r);
}

Complex RadialFlowSolver::sigma_rhs(double r, Complex sigma) const {
    Complex omega = wave_.omega;
    double f = bh_.f(r);
    double V = potential_V(r);

    // 系数
    Complex a = -1.0/r - V/(2.0*Complex(0, 1)*omega);
    Complex b = 2.0*Complex(0, 1)*omega/f + V/(Complex(0, 1)*omega);
    Complex c = -(-1.0/r + V/(2.0*Complex(0, 1)*omega));

    return a + b*sigma + c*sigma*sigma;
}

std::pair<double, Complex> RadialFlowSolver::sigma_initial(double epsilon) const {
    double r0 = bh_.r0;
    Complex omega = wave_.omega;
    int l = wave_.l;

    // 视界处：σ0 = 0
    Complex sigma0 = 0.0;

    // σ0' 的表达式
    double f0_prime = 1.0 / r0;
    Complex numerator = -f0_prime * (1.0 + static_cast<double>(l*(l+1))/(2.0*Complex(0,1)*omega*r0));
    Complex denominator = r0*f0_prime - 2.0*Complex(0,1)*omega*r0;
    Complex sigma0_prime = numerator / denominator;

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

SolveResult RadialFlowSolver::solve(double r_max, double epsilon, int n_points) const {
    auto [r_start, sigma_start] = sigma_initial(epsilon);

    double r_end = r_max * bh_.r0;
    double dr = (r_end - r_start) / static_cast<double>(n_points - 1);

    std::vector<double> r_vec;
    std::vector<Complex> sigma_vec;

    r_vec.reserve(n_points);
    sigma_vec.reserve(n_points);

    double r = r_start;
    Complex sigma = sigma_start;

    for (int i = 0; i < n_points; ++i) {
        r_vec.push_back(r);
        sigma_vec.push_back(sigma);

        if (i < n_points - 1) {
            rk4_step(r, sigma, dr);
            r += dr;
        }
    }

    // 反射系数 = -sigma (在无穷远处)
    Complex reflection_coeff = -sigma_vec.back();

    return {r_vec, sigma_vec, reflection_coeff, true};
}

Complex RadialFlowSolver::compute_reflection_coefficient(double r_max) const {
    SolveResult result = solve(r_max);
    return result.reflection_coeff;
}

}  // namespace radial_flow
