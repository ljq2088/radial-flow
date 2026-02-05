#!/usr/bin/env python3
"""
检查高频时的初值和 Riccati 方程系数
"""

import numpy as np
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)

print("=" * 70)
print("检查不同频率的初值和 Riccati 方程系数")
print("=" * 70)
print()

for omega in [0.5, 1.0, 2.0, 3.0, 5.0]:
    wave = WaveParams(omega=omega, l=2)
    solver = RadialFlowSolver(bh, wave)

    # 获取初值
    epsilon = 1e-6
    r0 = bh.r0

    # 计算 σ0'（从 .tex 公式）
    l = wave.l
    f0_prime = 1.0 / r0
    numerator = -f0_prime * (1.0 + (l*(l+1) + r0*f0_prime)/(2.0j*omega*r0))
    denominator = r0*f0_prime - 2.0j*omega*r0
    sigma0_prime = numerator / denominator

    r_start = r0 + epsilon * r0
    sigma_start = 0.0 + sigma0_prime * epsilon * r0

    print(f"ω = {omega}:")
    print(f"  σ0' = {sigma0_prime}")
    print(f"  |σ0'| = {np.abs(sigma0_prime):.6e}")
    print(f"  r_start = {r_start:.6f}")
    print(f"  σ_start = {sigma_start}")
    print(f"  |σ_start| = {np.abs(sigma_start):.6e}")

    # 计算 Riccati 方程在 r_start 处的系数
    r = r_start
    f = bh.f(r)
    f_prime = bh.f_prime(r)
    V = l*(l+1)/(r*r) + f_prime/r

    a = -1.0/r - V/(2.0j*omega)
    b = 2.0j*omega/f + V/(1.0j*omega)
    c = -(-1.0/r + V/(2.0j*omega))

    print(f"  Riccati 系数 (at r_start):")
    print(f"    a = {a}, |a| = {np.abs(a):.6e}")
    print(f"    b = {b}, |b| = {np.abs(b):.6e}")
    print(f"    c = {c}, |c| = {np.abs(c):.6e}")

    # 计算 dσ/dr
    dsigma_dr = a + b*sigma_start + c*sigma_start*sigma_start
    print(f"  dσ/dr = {dsigma_dr}")
    print(f"  |dσ/dr| = {np.abs(dsigma_dr):.6e}")
    print()

print("=" * 70)
