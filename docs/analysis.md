# 径向流方法分析与实现计划

## 项目架构

**C++ 核心 + Python 接口**：
- C++ 实现高性能计算核心
- pybind11 创建 Python 绑定
- Python 提供友好接口和可视化

## 推导验证

.tex 文件展示了 Schwarzschild 背景下标量场 (s=0) 的径向流方法。以下是我的验证：

### 1. 基本设置 (正确)

- Schwarzschild 度规：$f = 1 - r_0/r$，其中 $r_0 = 2M$
- 波动方程：$\partial_\mu(\sqrt{-g}g^{\mu\nu}\partial_\nu\psi) = 0$
- 分离变量：$\psi = R(r)Y_{lm}(\theta,\varphi)e^{-i\omega t}$

### 2. σ 的两种定义

文档给出了两种方法：

**第一种方法 (第31-50行)：**

$$\sigma = \frac{J}{-i\omega R}$$

其中 $J = g^{rr}\partial_r R$

- 边界条件：$\sigma_0 = 1$
- 这种定义在 R 的零点处有问题

**第二种方法 (第52-108行)：**

$$\sigma = \frac{J + i\omega R}{J - i\omega R}$$

这是**数值计算的首选定义**：
- 在无穷远处：$\sigma \sim -A/B$（反射系数）
- 视界边界条件：$\sigma_0 = 0$（纯入射波）

### 3. Riccati 方程 (第77行, 第99行)

最终的 ODE 为：

$$\partial_r\sigma = -\frac{1}{r} - \frac{V}{2i\omega} + \left[2i\omega f^{-1} + \frac{V}{i\omega}\right]\sigma - \left[-\frac{1}{r} + \frac{V}{2i\omega}\right]\sigma^2$$

其中势函数有两种形式：
- 完整形式：$V = r^{-2}l(l+1) + r^{-1}f'$（第77行）
- 简化形式：$V = r^{-2}l(l+1)$（第99行）

简化形式对应标准的 s=0 Regge-Wheeler 势。

### 4. 视界边界条件

在 $r = r_0$ 处：$f_0 = 0$，$f_0' = 1/r_0$

- $\sigma_0 = 0$
- $\sigma_0'$ 的表达式：

$$\sigma_0' = \frac{-f_0'\left(1 + \frac{l(l+1)}{2i\omega r_0}\right)}{r_0 f_0' - 2i\omega r_0} = \frac{-\frac{1}{r_0}\left(1 + \frac{l(l+1)}{2i\omega r_0}\right)}{1 - 2i\omega r_0}$$

这用于在视界附近启动积分。

---

## 实现计划

### 目录结构

```
radial_flow/
├── include/
│   └── radial_flow.hpp      # C++ 头文件
├── src/
│   └── radial_flow.cpp      # C++ 实现
├── python/
│   └── radial_flow.py       # Python 实现
├── tests/
│   └── test_scalar.py       # 验证测试
└── CMakeLists.txt           # 构建系统
```

### 核心组件

1. **ODE 系统**：实现 Riccati 方程 $\sigma' = f(r, \sigma)$
2. **边界条件**：在视界附近使用 Taylor 展开启动积分
3. **积分方法**：从 $r_0 + \epsilon$ 积分到大 $r$（使用 RK4 或自适应方法）
4. **输出**：从渐近 $\sigma$ 值提取 $-A/B$

### 数值注意事项

1. **视界奇点处理**：在 $r = r_0 + \epsilon$（如 $\epsilon = 10^{-6}r_0$）处开始积分
2. **Taylor 展开**：$\sigma(r_0 + \epsilon) \approx \sigma_0 + \sigma_0'\epsilon$
3. **积分终点**：选择足够大的 $r_{max}$（如 $100 r_0$）以达到渐近区域
4. **复数运算**：$\omega$ 可以是复数（用于准正规模）

---

## 当前进度

### ✅ 已完成

1. **C++ 核心库**
   - `include/radial_flow.hpp`：头文件（完整）
   - `src/radial_flow.cpp`：核心实现（完整）
   - `src/bindings.cpp`：pybind11 绑定（完整）

2. **Python 接口**
   - `python/radial_flow.py`：Python 包装层（完整）
   - `python/__init__.py`：包初始化（完整）

3. **构建系统**
   - `CMakeLists.txt`：CMake 配置（完整）

4. **文档**
   - `analysis.md`：推导验证和实现说明（本文件）

### 🔄 待完成

5. **构建脚本**
   - `setup.py`：Python 包安装脚本
   - `build.sh`：一键构建脚本

6. **测试和示例**
   - `tests/test_scalar.py`：单元测试
   - `examples/demo.py`：使用示例

7. **文档**
   - `README.md`：项目说明
   - `docs/usage.md`：使用指南

---

## 构建和使用步骤

### 1. 安装依赖

```bash
# 安装 pybind11
pip install pybind11

# 安装其他 Python 依赖
pip install numpy scipy matplotlib
```

### 2. 编译 C++ 核心

```bash
cd /home/ljq/code/radial_flow
mkdir build && cd build
cmake ..
make
make install  # 将 _radial_flow_cpp.so 安装到 python/ 目录
```

### 3. 运行 Python 程序

```bash
cd /home/ljq/code/radial_flow
python -m python.radial_flow
```

---

## 后续工作计划

### 第一阶段：完成基础功能（当前）
- [x] C++ 核心实现
- [x] Python 绑定
- [ ] 构建脚本
- [ ] 基础测试

### 第二阶段：验证和优化
- [ ] 与已知结果对比验证
- [ ] 性能优化（并行化、缓存等）
- [ ] 错误处理和边界情况

### 第三阶段：扩展到 s=-2
- [ ] 实现 Teukolsky 方程（s=-2）
- [ ] 引力波反射系数计算
- [ ] 更复杂的势函数

---

## 技术细节

### C++ 实现要点

1. **复数运算**：使用 `std::complex<double>`
2. **RK4 积分器**：4阶 Runge-Kutta 方法
3. **内存管理**：使用 `std::vector` 自动管理
4. **性能**：C++ 核心比纯 Python 快 10-100 倍

### Python 绑定要点

1. **pybind11**：自动类型转换
2. **numpy 兼容**：`std::vector` 自动转换为 numpy 数组
3. **复数支持**：`std::complex` 自动转换为 Python `complex`

### 数值稳定性

1. **视界附近**：使用 Taylor 展开避免奇点
2. **大 r 行为**：检查渐近收敛性
3. **步长选择**：自适应步长（未来可改进）
