# J 和 σ 的定义

## .tex 文件中的定义

### 第一种定义（第31行）❌ 不使用

```latex
σ := J / (-iωR)
```

其中 $J := g^{rr} \partial_r R$（第25行）

**问题**：这个定义在 R 的零点处有奇点。

---

### 第二种定义（第66-69行）✅ 使用这个

#### J 的定义（第66行）

```latex
J = (f/r) ∂_r(rR) = f ∂_r R + (f/r) R
```

**展开**：
$$J = f \frac{\partial R}{\partial r} + \frac{f}{r} R$$

其中 $f(r) = 1 - r_0/r$ 是度规函数。

**物理意义**：J 是径向流（radial flux）。

#### σ 的定义（第69行）

```latex
σ = (J + iωR) / (J - iωR)
```

**完整形式**：
$$\sigma = \frac{J + i\omega R}{J - i\omega R} = \frac{f\partial_r R + \frac{f}{r}R + i\omega R}{f\partial_r R + \frac{f}{r}R - i\omega R}$$

**渐近行为**（第69行）：
$$\sigma \sim -\frac{A_0}{B_0} \quad (r \to \infty)$$

---

## 在代码中的实现

### 当前代码中 J 的隐含定义

代码中没有显式定义 J，但从 Riccati 方程的推导可以看出，代码使用的是：

$$J = f \partial_r R + \frac{f}{r} R$$

这对应 .tex 第66行的定义。

### 当前代码中 σ 的定义

**代码位置**：`include/radial_flow.hpp` 第59-62行的注释

```cpp
/**
 * @brief 径向流方程（Riccati 方程）求解器
 *
 * 变量 sigma = (J + i*omega*R) / (J - i*omega*R) 满足 Riccati ODE。
 * 在无穷远处，sigma -> -A/B 给出反射系数。
 */
```

这对应 .tex 第69行的定义：
$$\sigma = \frac{J + i\omega R}{J - i\omega R}$$

---

## 关键性质

### 1. 视界处的边界条件

从 .tex 第83行：
$$\sigma_0 = 0$$

**物理意义**：视界处只有入射波，没有反射波。

**代码实现**：`src/radial_flow.cpp` 第36行
```cpp
Complex sigma0 = 0.0;
```

### 2. 无穷远处的渐近行为

从 .tex 第69行：
$$\sigma(r \to \infty) \sim -\frac{A_0}{B_0}$$

**物理意义**：σ 收敛到反射系数的负值。

**代码实现**：`src/radial_flow.cpp` 第84-85行
```cpp
// 反射系数 = -sigma (在无穷远处)
Complex reflection_coeff = -sigma_vec.back();
```

### 3. σ 的优点

相比第一种定义 $\sigma = J/(-i\omega R)$，第二种定义的优点：

1. **无奇点**：在 R 的零点处没有奇点
2. **明确的边界条件**：视界处 $\sigma_0 = 0$
3. **直接提取反射系数**：$A_0/B_0 = -\sigma(\infty)$

---

## Riccati 方程的推导

从 .tex 第72-77行，通过复杂的代数运算，可以从 σ 的定义推导出 Riccati 方程：

$$\frac{d\sigma}{dr} = -\frac{1}{r} - \frac{V}{2i\omega} + \left[\frac{2i\omega}{f} + \frac{V}{i\omega}\right]\sigma - \left[-\frac{1}{r} + \frac{V}{2i\omega}\right]\sigma^2$$

其中势函数（完整版本）：
$$V = \frac{l(l+1)}{r^2} + \frac{f'}{r}$$

**代码实现**：`src/radial_flow.cpp` 第17-28行

```cpp
Complex RadialFlowSolver::sigma_rhs(double r, Complex sigma) const {
    Complex omega = wave_.omega;
    double f = bh_.f(r);
    double V = potential_V(r);  // V = l(l+1)/r^2 + f'/r

    // 系数
    Complex a = -1.0/r - V/(2.0*Complex(0, 1)*omega);
    Complex b = 2.0*Complex(0, 1)*omega/f + V/(Complex(0, 1)*omega);
    Complex c = -(-1.0/r + V/(2.0*Complex(0, 1)*omega));

    return a + b*sigma + c*sigma*sigma;
}
```

---

## 总结

| 项目 | 定义 | 位置 |
|------|------|------|
| **J** | $J = f\partial_r R + \frac{f}{r}R$ | .tex 第66行 |
| **σ** | $\sigma = \frac{J + i\omega R}{J - i\omega R}$ | .tex 第69行 |
| **视界边界条件** | $\sigma_0 = 0$ | .tex 第83行 |
| **渐近行为** | $\sigma \sim -\frac{A_0}{B_0}$ | .tex 第69行 |
| **反射系数** | $\frac{A_0}{B_0} = -\sigma(\infty)$ | 代码第85行 |

**关键**：代码正确使用了 .tex 第66-69行的定义，这是数值稳定且物理清晰的定义。
