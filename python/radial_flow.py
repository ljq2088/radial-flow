"""
径向流计算：Schwarzschild 背景下的标量场 (s=0)

Python 接口层，调用 C++ 核心库进行计算。

参考：径向流(check) (1).tex
"""

import numpy as np
from typing import Optional
import matplotlib.pyplot as plt

# 导入 C++ 核心模块
try:
    from . import _radial_flow_cpp as cpp
except ImportError:
    import _radial_flow_cpp as cpp


class SchwarzschildParams:
    """Schwarzschild 黑洞参数（Python 包装）"""
    def __init__(self, r0: float = 2.0):
        self._cpp_obj = cpp.SchwarzschildParams()
        self._cpp_obj.r0 = r0

    @property
    def r0(self):
        return self._cpp_obj.r0

    @r0.setter
    def r0(self, value):
        self._cpp_obj.r0 = value

    def f(self, r):
        """度规函数 f = 1 - r0/r"""
        return self._cpp_obj.f(r)

    def f_prime(self, r):
        """f 的导数：f' = r0/r^2"""
        return self._cpp_obj.f_prime(r)


class WaveParams:
    """标量波参数（Python 包装）"""
    def __init__(self, omega: complex, l: int, m: int = 0):
        self._cpp_obj = cpp.WaveParams()
        self._cpp_obj.omega = omega
        self._cpp_obj.l = l
        self._cpp_obj.m = m

    @property
    def omega(self):
        return self._cpp_obj.omega

    @omega.setter
    def omega(self, value):
        self._cpp_obj.omega = value

    @property
    def l(self):
        return self._cpp_obj.l

    @l.setter
    def l(self, value):
        self._cpp_obj.l = value

    @property
    def m(self):
        return self._cpp_obj.m

    @m.setter
    def m(self, value):
        self._cpp_obj.m = value


class RadialFlowSolver:
    """
    径向流方程（Riccati 方程）求解器（Python 包装）

    调用 C++ 核心库进行高性能计算。
    """

    def __init__(self, bh: SchwarzschildParams, wave: WaveParams):
        self.bh = bh
        self.wave = wave
        # 创建 C++ 求解器对象
        self._cpp_solver = cpp.RadialFlowSolver(bh._cpp_obj, wave._cpp_obj)

    def solve(self, r_max: float = 200.0, epsilon: float = 1e-6,
              n_points: int = 10000) -> dict:
        """
        从视界积分到 r_max（调用 C++ 核心）

        参数：
            r_max: 积分终点（单位：r0）
            epsilon: 视界附近的起始偏移（相对于 r0）
            n_points: 输出点数

        返回：
            dict: 包含 r, sigma, reflection_coeff 等
        """
        # 调用 C++ 求解器
        cpp_result = self._cpp_solver.solve(r_max, epsilon, n_points)

        # 转换为 numpy 数组
        return {
            'r': np.array(cpp_result.r),
            'sigma': np.array(cpp_result.sigma),
            'reflection_coeff': cpp_result.reflection_coeff,
            'A_over_B': cpp_result.reflection_coeff,
            'success': cpp_result.success
        }

    def compute_reflection_coefficient(self, r_max: float = 200.0) -> complex:
        """计算反射系数 A/B（调用 C++ 核心）"""
        return self._cpp_solver.compute_reflection_coefficient(r_max)


def plot_sigma_evolution(solver: RadialFlowSolver, result: dict,
                         save_path: Optional[str] = None):
    """绘制 sigma 随 r 的演化"""
    r = result['r']
    sigma = result['sigma']
    r0 = solver.bh.r0

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # |sigma|
    axes[0, 0].plot(r/r0, np.abs(sigma))
    axes[0, 0].set_xlabel(r'$r/r_0$')
    axes[0, 0].set_ylabel(r'$|\sigma|$')
    axes[0, 0].set_title(r'$|\sigma|$ vs $r$')
    axes[0, 0].grid(True)

    # Re(sigma) 和 Im(sigma)
    axes[0, 1].plot(r/r0, sigma.real, label=r'Re($\sigma$)')
    axes[0, 1].plot(r/r0, sigma.imag, label=r'Im($\sigma$)')
    axes[0, 1].set_xlabel(r'$r/r_0$')
    axes[0, 1].set_ylabel(r'$\sigma$')
    axes[0, 1].set_title(r'$\sigma$ components vs $r$')
    axes[0, 1].legend()
    axes[0, 1].grid(True)

    # 复平面轨迹
    axes[1, 0].plot(sigma.real, sigma.imag)
    axes[1, 0].scatter([sigma.real[0]], [sigma.imag[0]], c='g', s=100,
                       label='Start (horizon)', zorder=5)
    axes[1, 0].scatter([sigma.real[-1]], [sigma.imag[-1]], c='r', s=100,
                       label='End (infinity)', zorder=5)
    axes[1, 0].set_xlabel(r'Re($\sigma$)')
    axes[1, 0].set_ylabel(r'Im($\sigma$)')
    axes[1, 0].set_title(r'$\sigma$ in complex plane')
    axes[1, 0].legend()
    axes[1, 0].grid(True)
    axes[1, 0].axis('equal')

    # 参数信息
    info_text = (f"Parameters:\n"
                 f"$\\omega$ = {solver.wave.omega}\n"
                 f"$l$ = {solver.wave.l}\n"
                 f"$r_0$ = {solver.bh.r0}\n\n"
                 f"Results:\n"
                 f"$A/B$ = {result['reflection_coeff']:.6f}\n"
                 f"$|A/B|$ = {np.abs(result['reflection_coeff']):.6f}")
    axes[1, 1].text(0.1, 0.5, info_text, transform=axes[1, 1].transAxes,
                    fontsize=12, verticalalignment='center',
                    family='monospace')
    axes[1, 1].axis('off')
    axes[1, 1].set_title('Summary')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"图像已保存到 {save_path}")

    plt.show()


def main():
    """主函数：演示径向流计算（使用 C++ 核心）"""
    # 设置参数
    bh = SchwarzschildParams(r0=2.0)  # M = 1
    wave = WaveParams(omega=0.5, l=2)

    print("=" * 60)
    print("径向流计算：Schwarzschild 标量场 (s=0)")
    print("使用 C++ 核心库进行高性能计算")
    print("=" * 60)
    print(f"黑洞参数：r0 = {bh.r0} (M = {bh.r0/2})")
    print(f"波参数：ω = {wave.omega}, l = {wave.l}")
    print()

    # 创建求解器并计算
    solver = RadialFlowSolver(bh, wave)
    result = solver.solve(r_max=200.0)

    print(f"积分成功：{result['success']}")
    print(f"反射系数 A/B = {result['reflection_coeff']}")
    print(f"|A/B| = {np.abs(result['reflection_coeff']):.6f}")
    print(f"arg(A/B) = {np.angle(result['reflection_coeff']):.6f} rad")
    print()

    # 绘图
    plot_sigma_evolution(solver, result)

    # 扫描不同频率
    print("\n频率扫描：")
    print("-" * 40)
    omegas = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0]
    for omega in omegas:
        wave_scan = WaveParams(omega=omega, l=2)
        solver_scan = RadialFlowSolver(bh, wave_scan)
        refl = solver_scan.compute_reflection_coefficient()
        print(f"ω = {omega:.1f}: A/B = {refl:.4f}, |A/B| = {np.abs(refl):.4f}")


if __name__ == "__main__":
    main()
