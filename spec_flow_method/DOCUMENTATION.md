# Schwarzschild 视界流与谱匹配方法

## 项目概述

本项目实现了求解 Schwarzschild 时空中标量场扰动的数值方法，采用**视界流（horizon flow）+ 谱方法匹配**的混合策略。

### 核心思想

1. **视界附近**：使用径向流方程（radial flow equation）从视界向外积分
2. **远场区域**：使用 Chebyshev 谱方法求解边值问题
3. **匹配点**：在中间位置 `z_match` 处连接两个解

---

## 物理背景

### 坐标系统

- 使用倒数坐标：`z = 1/r`
- 视界位置：`z_h = 1/r0`（其中 `r0` 是 Schwarzschild 半径）
- 无穷远：`z = 0`

### 控制方程

在 Schwarzschild 时空中，标量场扰动满足：

```
d²ψ/dr*² + ω²ψ = V(r)ψ
```

其中 `r*` 是龟坐标，`V(r)` 是有效势。

---

## 代码实现对比

### Mathematica Notebook vs Python 实现

#### 1. 径向流方程（Radial Flow）

**Mathematica 实现**（格林函数（齐次拼接）.nb）：
```mathematica
evz[σ_, z_] := {
  If[z == r0^(-1),
    (4ω²r0⁴ - l(l+1)r0²) / (-Iω - 2ω²r0),
    (-ω²z^(-2)/(1-r0*z) + l(l+1)) / (Iω*z²) - (Iω/(1-r0*z))σ²
  ],
  1
}
```
- 使用 Runge-Kutta 4 阶方法手动实现
- 从视界 `z = 1/r0` 开始，步长 `st = -0.00025`
- 迭代 1000 步

**Python 实现**（`flow.py`）：
```python
def sigma_flow_rhs(z: float, Sigma: complex, params: SolverParams) -> complex:
    r0 = params.r0
    w = params.omega
    ell = params.ell
    f = 1.0 - r0 * z
    term1 = (-w**2 * z**-2 / f + ell * (ell + 1.0)) / (1j * w * z**2)
    term2 = -(1j * w / f) * Sigma**2
    return term1 + term2
```
- 使用 `scipy.integrate.solve_ivp` 的 DOP853 方法（8阶 Runge-Kutta）
- 从 `z_h - flow_eps` 开始积分到 `z_match`
- 自适应步长，精度控制：`rtol=1e-10, atol=1e-12`

**一致性分析**：✅ **完全一致**
- 两者实现的是同一个 Riccati 方程
- Python 版本使用更高阶的积分器和自适应步长，精度更高

---

#### 2. Chebyshev 谱方法

**Mathematica 实现**：
- Notebook 中未明确展示 Chebyshev 网格构建代码
- 可能使用内置的谱方法函数

**Python 实现**（`chebyshev.py`）：
```python
def chebyshev_lobatto_grid(n: int, zc: float) -> ChebyshevGrid:
    k = np.arange(n, dtype=float)
    x = np.cos(np.pi * k / (n - 1))  # Gauss-Lobatto 节点

    # 构建微分矩阵
    c = np.ones(n, dtype=float)
    c[0] = 2.0
    c[-1] = 2.0
    c = c * ((-1.0) ** k)

    dX = x[:, None] - x[None, :]
    D_x = np.outer(c, 1.0 / c) / (dX + np.eye(n))
    D_x = D_x - np.diag(np.sum(D_x, axis=1))

    z = 0.5 * zc * (1.0 + x)  # 映射到 [0, zc]
    D_z = (2.0 / zc) * D_x

    return ChebyshevGrid(n=n, zc=zc, x=x, z=z, D=D_z)
```

**一致性分析**：✅ **标准实现**
- 使用标准的 Chebyshev-Gauss-Lobatto 节点
- 微分矩阵构建遵循经典算法
- 节点顺序：`z[0] = zc`（匹配点），`z[-1] = 0`（无穷远）

---

#### 3. 谱方法求解器

**Python 实现**（`spectral.py`）：

**正频模式算符**（`build_plus_operator`）：
```python
A = D @ diag(z²fz) @ D
A += (-2iω) * diag(1 - r0²z²) @ D
V = 2iωr0²z - r0z + (ω²r0²/fz)(2 - r0²z²) - ℓ(ℓ+1)
A += diag(V)
```

**负频模式算符**（`build_minus_operator`）：
```python
A = D @ diag(z²fz) @ D
A += (+2iω) * diag(1 - r0²z²) @ D
V = -2iωr0²z - r0z + (ω²r0²/fz)(2 - r0²z²) - ℓ(ℓ+1)
A += diag(V)
```

**边界条件**：
- 正频模式：`u_p(z=0) = 1`（无穷远归一化）
- 负频模式：在 `z=zc` 处与视界流解匹配

**一致性分析**：✅ **物理一致**
- 算符形式符合 Schwarzschild 时空的标量扰动方程
- 正负频模式通过 `±2iω` 项区分
- 边界条件正确实现了视界入射条件和无穷远渐近行为

---

#### 4. 匹配过程

**Python 实现**（`spectral.py`）：

```python
def boundary_plus_value(up, sigma_log, grid, params):
    zc = grid.zc
    d1u = grid.D[0, :] @ up  # 在 zc 处的导数
    u1 = up[0]
    val = (zc² * d1u + (zc - iω - iωr0*zc - zc²*sigma_log) * u1) * exp(iS(zc))
    return val

def solve_um(grid, sigma_log, bdry_plus, params):
    A = build_minus_operator(grid, params)
    rhs = np.zeros(grid.n, dtype=complex)
    rhs[0] = bdry_plus
    A[0, :] = boundary_minus_row(sigma_log, grid, params)
    return np.linalg.solve(A, rhs)
```

**匹配逻辑**：
1. 从视界流得到 `σ_log = d(log ψ)/dz` 在 `z_match` 处的值
2. 求解正频模式 `u_p`，在无穷远归一化
3. 计算 `u_p` 在匹配点的"边界值"（包含相位因子）
4. 用此边界值作为负频模式 `u_m` 的边界条件
5. 计算振幅比 `ratio = u_p(z=0) / u_m(z=0)`

**一致性分析**：✅ **方法正确**
- 匹配条件正确连接了视界流解和谱方法解
- 相位因子 `exp(iS)` 正确处理了 WKB 相位
- 振幅比给出了反射/透射系数

---


## 核心算法详解

### 1. 视界流方程

从视界附近开始，定义辅助变量 `Σ = (1-r0*z)/(iω) * dψ/dz / ψ`，满足 Riccati 方程：

```
dΣ/dz = [(-ω²/z²)/(1-r0*z) + ℓ(ℓ+1)] / (iω*z²) - (iω/(1-r0*z)) * Σ²
```

**视界边界条件**：
```python
Σ'(z_h) = [4ω²r0⁴ - ℓ(ℓ+1)r0²] / [-iω - 2ω²r0]
```

**实现细节**：
- 从 `z_start = z_h - ε` 开始（避免奇点）
- 线性外推：`Σ(z_start) = r0² + Σ'(z_h) * (z_start - z_h)`
- 积分到匹配点 `z_match`
- 转换为对数导数：`σ_log = iω*Σ/(1-r0*z)`

### 2. Chebyshev 谱方法

**微分方程**：
```
d/dz[z²f(z) dψ/dz] + [V(z) ± 2iω(1-r0²z²) d/dz]ψ = 0
```

其中 `f(z) = 1 - r0*z`，`±` 对应正/负频模式。

**离散化**：
- 使用 Chebyshev 微分矩阵 `D` 近似导数
- 将微分方程转化为线性系统 `A·u = rhs`
- 边界条件通过修改矩阵第一行实现

**正频模式边界条件**：
- `u_p(z=0) = 1`（无穷远归一化）

**负频模式边界条件**：
- 在 `z=zc` 处匹配视界流解
- 匹配条件：`[z²du/dz + (z - iω - iωr0*z - z²*σ_log)u] * exp(iS) = bdry_plus`

### 3. 相位函数

```python
S(z) = ω * (1/z + r0 * log(1/(r0*z)))
```

这是 WKB 近似中的相位积分，用于连接不同区域的解。


## 代码结构

```
src/schw_matcher/
├── __init__.py          # 包导出
├── params.py            # 参数定义（SolverParams）
├── flow.py              # 视界流积分
├── chebyshev.py         # Chebyshev 网格和微分矩阵
├── spectral.py          # 谱方法求解器
├── matching.py          # 主匹配流程
└── solver.py            # 高层接口

examples/
└── run_case.py          # 使用示例

tests/
└── smoke_test.py        # 基本测试
```

## 使用方法

### 基本用法

```python
from schw_matcher import SolverParams, solve_case

# 设置参数
params = SolverParams(
    r0=2.0,           # Schwarzschild 半径
    omega=0.2,        # 频率
    ell=0,            # 角量子数
    z_match=0.25,     # 匹配点位置
    n_cheb=25,        # Chebyshev 节点数
    flow_eps=1e-6     # 视界偏移量
)

# 求解
result = solve_case(params)

# 提取结果
ratio = result.spectral.ratio_up_over_um
amplitude = abs(ratio)
phase = result.amplitude_ratio_arg

print(f"振幅比: {amplitude}")
print(f"相位: {phase}")
```

### 参数说明

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `r0` | float | 2.0 | Schwarzschild 半径 |
| `omega` | float | 0.2 | 频率（必须非零） |
| `ell` | int | 0 | 角量子数 ℓ |
| `z_match` | float | 0.25 | 匹配点（需满足 0 < z_match < 1/r0） |
| `n_cheb` | int | 25 | Chebyshev 节点数（≥4） |
| `flow_eps` | float | 1e-6 | 视界偏移量 |
| `flow_rtol` | float | 1e-10 | 流方程相对误差 |
| `flow_atol` | float | 1e-12 | 流方程绝对误差 |


## 结果解释

### MatchingResult 对象

```python
@dataclass
class MatchingResult:
    params: SolverParams              # 输入参数
    flow: FlowResult                  # 视界流结果
    grid: ChebyshevGrid               # Chebyshev 网格
    spectral: SpectralState           # 谱方法结果

    @property
    def amplitude_ratio_abs(self) -> float:
        """振幅比的模"""

    @property
    def amplitude_ratio_arg(self) -> float:
        """振幅比的幅角"""
```

### FlowResult 对象

```python
@dataclass
class FlowResult:
    z_start: float                    # 起始位置
    z_end: float                      # 终止位置
    Sigma_start: complex              # 起始 Σ 值
    Sigma_end: complex                # 终止 Σ 值
    sigma_logderiv: complex           # 对数导数 σ
    nfev: int                         # 函数评估次数
    success: bool                     # 是否成功
    message: str                      # 状态信息
```

### SpectralState 对象

```python
@dataclass
class SpectralState:
    grid: ChebyshevGrid               # 网格
    fz: ArrayC                        # f(z) = 1 - r0*z
    up: ArrayC                        # 正频模式解
    um: ArrayC                        # 负频模式解
    bdry_plus: complex                # 边界值
    ratio_up_over_um: complex         # 振幅比 u_p(0)/u_m(0)
```

