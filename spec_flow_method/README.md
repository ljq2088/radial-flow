# Schwarzschild horizon-flow + spectral matching

这是把 Mathematica notebook `格林函数（齐次拼接）.nb` 里的“视界到匹配点 Riccati flow + 匹配点到无穷远 Chebyshev 拼接”整理成 Python 的最小可运行版本。

## 数学结构

令

\[
z = 1/r,\qquad f(z)=1-r_0 z,\qquad z_h = 1/r_0.
\]

原始径向函数 \(R(z)\) 满足

\[
\partial_z\!\bigl(f\,\partial_z R\bigr)
+\left(
\frac{\omega^2}{f z^4}-\frac{l(l+1)}{z^2}
\right)R=0.
\]

内侧 flow 用的变量是

\[
\Sigma(z)=\frac{1-r_0 z}{i\omega}\frac{R'(z)}{R(z)},
\]

它满足 Riccati 方程

\[
\Sigma'
=
\frac{-\omega^2 z^{-2}/(1-r_0 z)+l(l+1)}{i\omega z^2}
-
\frac{i\omega}{1-r_0 z}\Sigma^2.
\]

在视界纯入射边界下，\(\Sigma(z_h)=r_0^2\)。推进到匹配点 \(z_1\) 后，真实对数导数是

\[
\sigma_1 = \left.\frac{R'}{R}\right|_{z_1}
=\frac{i\omega}{1-r_0 z_1}\Sigma(z_1).
\]

外侧谱方法改用 \(\Phi=rR\)，即 \(R=z\Phi\)。再把 \(\Phi\) 写成两组 Jost 型基解

\[
\Phi_\pm(z)=e^{\pm iS(z)}u_\pm(z),
\qquad
S(z)=\omega\left(\frac1z+r_0\log\frac1{r_0 z}\right).
\]

最后在 \(z_1\) 处施加

\[
R'(z_1)=\sigma_1 R(z_1)
\qquad \Longleftrightarrow \qquad
(z\Phi)'_{z_1}=\sigma_1 z_1\Phi(z_1),
\]

并用这个条件把 \(u_+\) 和 \(u_-\) 拼起来。

## 目录

```text
schw_matcher_project/
├── README.md
├── pyproject.toml
├── requirements.txt
├── src/
│   └── schw_matcher/
│       ├── __init__.py
│       ├── params.py
│       ├── flow.py
│       ├── chebyshev.py
│       ├── spectral.py
│       ├── matching.py
│       └── solver.py
├── examples/
│   └── run_case.py
└── tests/
    └── smoke_test.py
```

## 运行

```bash
pip install -r requirements.txt
export PYTHONPATH=./src
python examples/run_case.py
```

## 当前定位

这是一份“结构清楚、和 notebook 对应”的最小原型。下一步最值得做的是：

1. 视界初值提升到更高阶近视界展开；
2. 对谱离散的条件数与收敛性做系统检查；
3. 把 \(u_+\)、\(u_-\) 的归一化和最终 Jost 系数的定义再规范化。
