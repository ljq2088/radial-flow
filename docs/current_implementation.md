# 当前代码实现的公式分析

## 当前代码使用的公式

### 1. σ 的定义

```cpp
// 代码中隐含的定义
σ = (J + iωR) / (J - iωR)
```

其中：
- $R(r)$ 是径向波函数
- $J = f(r) \partial_r R$ 是径向流
- $\omega$ 是频率
- $f(r) = 1 - r_0/r$ 是度规函数

### 2. 势函数 V（当前实现）

**代码位置**：`src/radial_flow.cpp` 第12-15行

```cpp
double RadialFlowSolver::potential_V(double r) const {
    int l = wave_.l;
    return static_cast<double>(l * (l + 1)) / (r * r);
}
```

**对应公式**：
$$V(r) = \frac{l(l+1)}{r^2}$$

**来源**：.tex 文件第99行（简化版本）

### 3. Riccati 方程（当前实现）

**代码位置**：`src/radial_flow.cpp` 第17-28行

```cpp
Complex RadialFlowSolver::sigma_rhs(double r, Complex sigma) const {
    Complex omega = wave_.omega;
    double f = bh_.f(r);
    double V = potential_V(r);  // V = l(l+1)/r^2

    // 系数
    Complex a = -1.0/r - V/(2.0*Complex(0, 1)*omega);
    Complex b = 2.0*Complex(0, 1)*omega/f + V/(Complex(0, 1)*omega);
    Complex c = -(-1.0/r + V/(2.0*Complex(0, 1)*omega));

    return a + b*sigma + c*sigma*sigma;
}
```

**对应公式**：
$$\frac{d\sigma}{dr} = a + b\sigma + c\sigma^2$$

其中：
$$a = -\frac{1}{r} - \frac{V}{2i\omega}$$

$$b = \frac{2i\omega}{f} + \frac{V}{i\omega}$$

$$c = \frac{1}{r} - \frac{V}{2i\omega}$$

代入 $V = \frac{l(l+1)}{r^2}$：

$$\frac{d\sigma}{dr} = -\frac{1}{r} - \frac{l(l+1)}{2i\omega r^2} + \left[\frac{2i\omega}{f} + \frac{l(l+1)}{i\omega r^2}\right]\sigma - \left[-\frac{1}{r} + \frac{l(l+1)}{2i\omega r^2}\right]\sigma^2$$

**来源**：.tex 文件第99行（简化版本）

### 4. 视界边界条件（当前实现）

**代码位置**：`src/radial_flow.cpp` 第30-48行

```cpp
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
```

**对应公式**：
$$\sigma_0 = 0$$

$$\sigma_0' = \frac{-f_0' \left(1 + \frac{l(l+1)}{2i\omega r_0}\right)}{r_0 f_0' - 2i\omega r_0}$$

其中 $f_0' = \frac{1}{r_0}$，代入后：

$$\sigma_0' = \frac{-\frac{1}{r_0} \left(1 + \frac{l(l+1)}{2i\omega r_0}\right)}{1 - 2i\omega r_0}$$

**来源**：.tex 文件第105行（简化版本）

### 5. 渐近行为（当前实现隐含）

**代码位置**：`src/radial_flow.cpp` 第84-85行

```cpp
// 反射系数 = -sigma (在无穷远处)
Complex reflection_coeff = -sigma_vec.back();
```

这意味着代码假设：
$$\sigma(r \to \infty) \to -\frac{A}{B}$$

**对应的渐近行为**：
$$R(r) \sim A e^{i\omega r} + B e^{-i\omega r} \quad (r \to \infty)$$

**来源**：.tex 文件第89行（简化版本）

---

## 问题诊断

### 观察到的现象

从图像 `sigma_complex_plane_omega1.0.png` 可以看到：
- **σ 在复平面上做圆周运动（震荡）**
- **σ 不收敛到固定点**
- 这导致反射系数的提取不准确

### 问题根源

当前代码使用的是 **.tex 文件中的简化版本**（第87-108行），对应的渐近行为：
$$R \sim A e^{i\omega r} + B e^{-i\omega r}$$

这个渐近行为**不正确**，因为它忽略了：
1. **龟坐标** $r_*$ 的影响
2. **$1/r$ 衰减因子**
3. **$f'$ 项在势函数中的贡献**

### 正确的渐近行为

根据 .tex 文件第52-86行（完整版本），正确的渐近行为应该是：
$$R \sim \left(A_0 + \frac{A_1}{r}\right) r^{-1} e^{i\omega r_*} + \left(B_0 + \frac{B_1}{r}\right) r^{-1} e^{-i\omega r_*}$$

其中龟坐标：
$$r_* = r + r_0 \ln\left(\frac{r}{r_0} - 1\right)$$

---

## 需要修改的地方

### 修改 1：势函数

**当前**：
```cpp
return static_cast<double>(l * (l + 1)) / (r * r);
```
$$V = \frac{l(l+1)}{r^2}$$

**应该改为**：
```cpp
double f_prime = bh_.f_prime(r);
return static_cast<double>(l * (l + 1)) / (r * r) + f_prime / r;
```
$$V = \frac{l(l+1)}{r^2} + \frac{f'}{r}$$

**理由**：完整版本的 Riccati 方程（.tex 第77行）包含 $r^{-1}f'$ 项

### 修改 2：边界条件

**当前**：
```cpp
Complex numerator = -f0_prime * (1.0 + static_cast<double>(l*(l+1))/(2.0*Complex(0,1)*omega*r0));
```

**应该改为**：
```cpp
Complex numerator = -f0_prime * (1.0 + static_cast<double>(l*(l+1) + r0*f0_prime)/(2.0*Complex(0,1)*omega*r0));
```

**理由**：完整版本的边界条件（.tex 第83行）分子中是 $l(l+1) + r_0 f_0'$

---

## 公式对比表

| 项目 | 当前实现（简化版本） | 应该使用（完整版本） | .tex 位置 |
|------|---------------------|---------------------|-----------|
| 渐近行为 | $R \sim A e^{i\omega r} + B e^{-i\omega r}$ | $R \sim (A_0 + A_1/r) r^{-1} e^{i\omega r_*} + (B_0 + B_1/r) r^{-1} e^{-i\omega r_*}$ | 第55行 vs 第89行 |
| 势函数 $V$ | $\frac{l(l+1)}{r^2}$ | $\frac{l(l+1)}{r^2} + \frac{f'}{r}$ | 第99行 vs 第77行 |
| 边界条件分子 | $l(l+1)$ | $l(l+1) + r_0 f_0'$ | 第105行 vs 第83行 |
| Riccati 方程 | 不含 $f'/r$ 项 | 含 $f'/r$ 项 | 第99行 vs 第77行 |

---

## 预期修复效果

修复后，应该观察到：
1. **σ 平滑收敛**：在复平面上从原点平滑收敛到固定点 $-A_0/B_0$，无震荡
2. **物理合理的反射系数**：
   - 低频：$|A/B| \to 1$（完全反射）
   - 高频：$|A/B| \to 0$（完全透射）
   - 中间频率：平滑单调过渡
3. **无圆周运动**：σ 的轨迹应该是平滑曲线，不是圆

---

## 参考文献

- **径向流(check) (1).tex**：
  - 第52-86行：完整版本（正确）✅
  - 第87-108行：简化版本（当前使用，错误）❌
