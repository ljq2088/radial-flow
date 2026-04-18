# Kerr Teukolsky Riccati + Spectral Matching

这是一个仿照之前 Schwarzschild 项目搭建的 Python 项目，用于实现：

1. Kerr 背景下 `s=-2` Teukolsky 径向方程的内侧 Riccati flow；
2. 从有限匹配点 `r_match` 到无穷远的外侧 Chebyshev 谱方法；
3. 在匹配点处用 `(R_z / R)|_{z_1}` 进行拼接。

## 当前实现采用的定义

内侧变量：

\[
J = \frac{\Delta}{r^2+a^2} R_r,
\qquad
\sigma = \frac{J}{i k_H R},
\qquad
k_H = \omega - m \Omega_H.
\]

对应 Riccati 方程：

\[
\frac{d\sigma}{dr}
=
\left(2\frac{\Delta'}{\Delta} - \frac{2r}{r^2+a^2}\right)\sigma
- i k_H \frac{r^2+a^2}{\Delta}\sigma^2
- \frac{V(r)}{i k_H (r^2+a^2)}.
\]

外侧谱方法在 `z = 1/r` 上工作，区间为 `z in [0, z_match]`，并使用

\[
R_+(z)=A_+(z)u_+(z), \qquad A_+(z)=z^{-3} e^{+i\omega r_*(z)},
\]
\[
R_-(z)=A_-(z)u_-(z), \qquad A_-(z)=z e^{-i\omega r_*(z)}.
\]

## 目录结构

```text
kerr_matcher_project/
├── README.md
├── pyproject.toml
├── requirements.txt
├── src/
│   └── kerr_matcher/
│       ├── __init__.py
│       ├── params.py
│       ├── kerr.py
│       ├── chebyshev.py
│       ├── flow.py
│       ├── spectral.py
│       ├── matching.py
│       └── solver.py
├── examples/
│   └── run_case.py
└── tests/
    └── smoke_test.py
```

## 安装与运行

在项目根目录下：

```bash
pip install -r requirements.txt
PYTHONPATH=./src python examples/run_case.py
PYTHONPATH=./src python tests/smoke_test.py
```

或安装成 editable package：

```bash
pip install -e .
python examples/run_case.py
```

## 关于 separation constant

Kerr 下真正的角向 separation constant 一般应来自自旋加权球/椭球谐函数本征值。
当前实现里：

- 若提供 `lambda_sep`，则直接使用；
- 若不提供，则默认使用 Schwarzschild 极限近似
  \(
  \lambda = \ell(\ell+1)-2
  \)
  作为占位。

因此，这个项目当前的定位是：

- 结构完整；
- 方便你继续替换更精确的 `lambda_sep`；
- 方便你继续改进近视界启动、外侧算子和收敛性测试。

## 当前数值实现的诚实说明

这一版是“可运行骨架 + 与你当前推导一致的模块化实现”，重点是：

- 保持与推导中 `J, sigma, R_\pm = A_\pm u_\pm` 的记号一致；
- 外侧算子使用已正则化后的有限系数 `B_\pm(z), C_\pm(z)`；
- 便于后续继续做：
  - 更高阶近视界展开；
  - 真正的 spheroidal eigenvalue；
  - 不同 `r_match / n_cheb / flow_eps` 的收敛性扫描。
