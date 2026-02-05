# 径向流项目完整文档

## 目录

1. [项目概述](#项目概述)
2. [环境设置](#环境设置)
3. [构建过程](#构建过程)
4. [计算方法](#计算方法)
5. [使用指南](#使用指南)
6. [测试结果](#测试结果)
7. [项目结构](#项目结构)
8. [故障排除](#故障排除)
9. [参考文献](#参考文献)

---

## 项目概述

### 项目目标

本项目实现了**径向流方法**（Radial Flow Method）用于计算 Schwarzschild 黑洞背景下标量场扰动的**反射系数**。这是黑洞扰动理论中的一个基础问题，对理解黑洞的散射性质和准正规模具有重要意义。

### 物理背景

#### Schwarzschild 黑洞

Schwarzschild 黑洞是最简单的黑洞解，其度规为：

$$ds^2 = -f(r)dt^2 + f(r)^{-1}dr^2 + r^2(d\theta^2 + \sin^2\theta d\varphi^2)$$

其中度规函数：

$$f(r) = 1 - \frac{r_0}{r} = 1 - \frac{2M}{r}$$

- $r_0 = 2M$ 是 Schwarzschild 半径（视界位置）
- $M$ 是黑洞质量（本项目中取 $M = 1$，即 $r_0 = 2$）

#### 标量场扰动（s=0）

考虑标量场 $\psi$ 在 Schwarzschild 背景下满足的波动方程：

$$\nabla_\mu\nabla^\mu\psi = 0$$

通过分离变量：

$$\psi(t, r, \theta, \varphi) = e^{-i\omega t} R(r) Y_{lm}(\theta, \varphi)$$

可以将问题简化为径向方程。这里：
- $\omega$ 是频率
- $l, m$ 是角动量量子数
- $Y_{lm}$ 是球谐函数

#### 反射系数的物理意义

当波从无穷远入射到黑洞时，部分被黑洞吸收，部分被反射回无穷远。在无穷远处，径向波函数的渐近行为为：

$$R(r) \sim B e^{i\omega r_*} + A e^{-i\omega r_*} \quad (r \to \infty)$$

其中：
- $B$ 是入射波振幅
- $A$ 是反射波振幅
- $r_* = r + r_0 \ln|r/r_0 - 1|$ 是龟坐标

**反射系数**定义为：

$$\mathcal{R} = \frac{A}{B}$$

其模 $|\mathcal{R}|$ 表示反射波与入射波的振幅比：
- $|\mathcal{R}| = 1$：完全反射（低频极限）
- $|\mathcal{R}| = 0$：完全透射（高频极限）

### 技术架构

本项目采用 **C++ 核心 + Python 接口** 的混合架构：

- **C++ 核心库**：实现高性能数值计算（RK4 积分器）
- **pybind11 绑定**：在 C++ 和 Python 之间建立桥梁
- **Python 接口层**：提供友好的 API 和可视化功能

这种架构结合了 C++ 的计算性能和 Python 的易用性。

---

## 环境设置

### 系统要求

- **操作系统**：Linux 或 WSL2（Windows Subsystem for Linux 2）
- **编译器**：g++ 11.4.0 或更高版本（支持 C++17）
- **CMake**：3.15 或更高版本
- **Python**：3.8 或更高版本

### 依赖列表

#### Python 依赖

| 包名 | 版本 | 用途 |
|------|------|------|
| pybind11 | 3.0.1 | C++/Python 绑定 |
| numpy | 2.3.3 | 数值计算 |
| scipy | 1.16.2 | 科学计算（未来扩展） |
| matplotlib | 3.10.6 | 可视化 |

#### 系统依赖

- **CMake** 3.22.1：构建系统
- **g++** 11.4.0：C++ 编译器
- **Python 开发头文件**：python3-dev

### 安装步骤

#### 1. 安装系统依赖

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install cmake g++ python3-dev

# 验证安装
cmake --version
g++ --version
python3 --version
```

#### 2. 安装 Python 依赖

```bash
# 使用 pip 安装
pip install pybind11 numpy scipy matplotlib

# 或使用 requirements.txt（如果提供）
pip install -r requirements.txt
```

#### 3. 验证 pybind11 安装

```bash
python3 -c "import pybind11; print(pybind11.get_cmake_dir())"
```

应该输出 pybind11 的 CMake 配置目录路径。

---

## 构建过程

### 快速构建（推荐）

项目提供了一键构建脚本：

```bash
cd /home/ljq/code/radial_flow
./build.sh
```

### 手动构建步骤

如果需要更多控制，可以手动执行构建步骤：

#### 1. 创建构建目录

```bash
cd /home/ljq/code/radial_flow
mkdir -p build
cd build
```

#### 2. 配置 CMake

```bash
cmake ..
```

CMake 会：
- 查找 pybind11
- 配置 C++17 标准
- 设置包含目录
- 配置构建目标

#### 3. 编译

```bash
make -j$(nproc)
```

这会编译：
- `libradial_flow_core.a`：C++ 静态库
- `_radial_flow_cpp.so`：Python 扩展模块

#### 4. 安装

```bash
make install
```

这会将 `_radial_flow_cpp.so` 安装到 `python/` 目录。

### 构建输出说明

成功构建后，你应该看到：

```
radial_flow/
├── build/
│   ├── libradial_flow_core.a      # C++ 静态库
│   └── _radial_flow_cpp.so        # Python 扩展模块（临时）
└── python/
    └── _radial_flow_cpp.so        # Python 扩展模块（安装后）
```

### 验证构建

```bash
cd /home/ljq/code/radial_flow
python3 -c "from python import _radial_flow_cpp; print('Build successful!')"
```

如果没有错误，说明构建成功。

---

## 计算方法

### 理论基础

#### 径向流方法概述

径向流方法通过引入辅助变量 $\sigma$ 将二阶径向方程转化为一阶 Riccati 方程，然后从视界向外积分。

#### Riccati 变量定义

定义：

$$\sigma = \frac{J + i\omega R}{J - i\omega R}$$

其中：
- $R(r)$ 是径向波函数
- $J = g^{rr} \partial_r R = f(r) R'(r)$ 是径向流
- $\omega$ 是频率

这个定义的优点：
1. 在视界处有明确的边界条件：$\sigma_0 = 0$
2. 在无穷远处：$\sigma \to -A/B$（反射系数）
3. 避免了 $R$ 的零点问题

