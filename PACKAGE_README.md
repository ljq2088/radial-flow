# Radial Flow 计算包

Schwarzschild 黑洞背景下标量场散射的径向流方法计算包

## 📋 项目简介

本项目实现了使用径向流方法计算 Schwarzschild 黑洞背景下标量场的反射系数 A/B。通过从视界向外积分 Riccati 方程，可以高效地计算散射振幅。

### 主要特性

- ✅ **高性能 C++ 核心**：使用 RK4 方法进行数值积分
- ✅ **Python 接口**：通过 pybind11 提供易用的 Python API
- ✅ **两种势函数版本**：
  - 简化版本：V = l(l+1)/r²
  - 完整版本：V = l(l+1)/r² + f'/r（包含次领头项）
- ✅ **自适应收敛检测**：基于相位修正的智能收敛判断
- ✅ **可视化支持**：生成详细的结果图表

## 📁 文件结构

```
radial_flow/
├── include/
│   └── radial_flow.hpp          # C++ 头文件
├── src/
│   ├── radial_flow.cpp          # C++ 实现
│   └── bindings.cpp             # Python 绑定
├── python/
│   └── radial_flow.py           # Python 包装层
├── build/                       # 编译输出目录
├── figures/                     # 图表输出目录
├── CMakeLists.txt               # CMake 配置文件
├── demo.py                      # 主演示脚本
├── quick_test.py                # 快速测试脚本
├── example_usage.py             # 使用示例
└── PACKAGE_README.md            # 本文件
```

## 🔧 依赖项

### 系统要求
- Linux / macOS / WSL
- C++17 或更高版本
- CMake 3.15+
- Python 3.8+

### C++ 依赖
- pybind11 (自动下载)
- 标准 C++ 库

### Python 依赖
```bash
pip install numpy matplotlib
```

## 🚀 编译和安装

### 步骤 1: 克隆或下载项目

```bash
cd /path/to/radial_flow
```

### 步骤 2: 创建构建目录

```bash
mkdir -p build
cd build
```

### 步骤 3: 配置和编译

```bash
cmake ..
make
```

### 步骤 4: 复制编译产物到 Python 目录

```bash
cp _radial_flow_cpp.*.so ../python/
```

或者使用一键编译脚本：

```bash
./build.sh
```

## 📖 使用方法

### 方法 1: 运行演示脚本（推荐）

```bash
python3 demo.py
```

这将：
1. 计算简化版本和完整版本的反射系数
2. 输出详细的文本结果
3. 生成对比可视化图表（保存到 `figures/radial_flow_comparison.png`）

### 方法 2: 快速测试

```bash
python3 quick_test.py
```

快速验证安装是否成功，并查看基本结果。

### 方法 3: 在 Python 代码中调用

```python
import sys
sys.path.insert(0, '/path/to/radial_flow')

from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

# 1. 设置黑洞参数
bh = SchwarzschildParams(r0=2.0)  # Schwarzschild 半径

# 2. 设置波参数
wave = WaveParams(
    omega=0.2,                    # 频率
    l=0,                          # 角动量量子数
    use_subleading_term=True      # True=完整版本, False=简化版本
)

# 3. 创建求解器并计算
solver = RadialFlowSolver(bh, wave)
result = solver.solve(r_max=200.0)

# 4. 获取结果
A_over_B = result['reflection_coeff']
print(f"|A/B| = {abs(A_over_B):.6f}")
```

## 🎯 核心 API

### SchwarzschildParams

黑洞参数类。

```python
bh = SchwarzschildParams(r0=2.0)  # r0 = 2M
```

**属性：**
- `r0`: Schwarzschild 半径（默认 2.0）

**方法：**
- `f(r)`: 度规函数 f = 1 - r₀/r
- `f_prime(r)`: f 的导数 f' = r₀/r²

### WaveParams

波参数类。

```python
wave = WaveParams(omega=0.2, l=0, use_subleading_term=True)
```

**参数：**
- `omega` (complex): 频率（可以是复数，用于准正规模）
- `l` (int): 角动量量子数（l ≥ 0）
- `m` (int): 方位角量子数（默认 0，Schwarzschild 不依赖 m）
- `use_subleading_term` (bool):
  - `True`: 完整版本，V = l(l+1)/r² + f'/r
  - `False`: 简化版本，V = l(l+1)/r²

### RadialFlowSolver

径向流求解器。

```python
solver = RadialFlowSolver(bh, wave)
result = solver.solve(r_max=200.0, epsilon=1e-6, n_points=10000)
```

**solve() 参数：**
- `r_max` (float): 最小积分终点（单位：r₀，默认 200.0）
- `epsilon` (float): 视界偏移参数（默认 1e-6）
- `n_points` (int): 输出点数（默认 10000）
- `convergence_threshold` (float): 收敛阈值（默认 0.01 即 1%）
- `n_periods` (int): 检测周期数（默认 2）
- `max_r_factor` (float): 最大扩展倍数（默认 5.0）

**返回值：**
- `r`: 半径数组
- `sigma`: 修正后的 σ 数组
- `reflection_coeff`: 反射系数 A/B
- `success`: 是否成功收敛

## 📊 示例输出

### 文本输出

```
================================================================================
Schwarzschild 黑洞标量场散射计算
================================================================================
计算时间: 2026-02-04 19:30:00

黑洞参数:
  Schwarzschild 半径: r₀ = 2.0
  质量: M = 1.0

波参数:
  频率: ω = 0.2
  角动量量子数: l = 0

--------------------------------------------------------------------------------

[1/2] 计算简化版本 (V = l(l+1)/r²)...
[2/2] 计算完整版本 (V = l(l+1)/r² + f'/r)...

================================================================================
计算结果
================================================================================

【简化版本】V = l(l+1)/r²
  势函数: 仅包含角动量势垒
  收敛状态: False
  积分终点: r = 800.03
  反射系数: A/B = -0.263929 + -0.151991i
  模: |A/B| = 0.304564
  相位: arg(A/B) = -2.619099 rad

【完整版本】V = l(l+1)/r² + f'/r
  势函数: 包含角动量势垒 + 时空曲率修正
  收敛状态: True
  积分终点: r = 400.00
  反射系数: A/B = -0.213005 + -0.099159i
  模: |A/B| = 0.234955
  相位: arg(A/B) = -2.705903 rad

【对比分析】
  |A/B| 差异: 0.069610
  相对差异: 22.86%
  物理解释: 次领头项 f'/r 使势垒增高，导致反射系数减小约 23%
```

### 可视化输出

运行 `demo.py` 会生成包含以下内容的图表：
1. 简化版本的 |σ| vs r
2. 完整版本的 |σ| vs r
3. 两个版本的对比
4. 简化版本的复平面轨迹
5. 完整版本的复平面轨迹
6. 结果总结面板

图表保存在 `figures/radial_flow_comparison.png`

## 🔬 物理背景

### 径向流方法

径向流方法通过引入变量 σ = (J + iωR) / (J - iωR) 将二阶 ODE 转化为一阶 Riccati 方程：

```
dσ/dr = a + bσ + cσ²
```

其中系数 a, b, c 取决于势函数 V(r) 和度规函数 f(r)。

### 两种势函数版本

**简化版本：** V = l(l+1)/r²
- 仅考虑角动量势垒
- 对应简化的渐近行为

**完整版本：** V = l(l+1)/r² + f'/r
- 包含时空曲率的次领头修正
- 对应正确的渐近行为：R ~ (A₀ + A₁/r) r⁻¹ e^(iωr*)

### 边界条件

- **视界处 (r = r₀)**：纯入射波条件 σ₀ = 0
- **无穷远处 (r → ∞)**：σ → -A/B（反射系数）

## 🐛 故障排除

### 编译错误

**问题：** `pybind11 not found`
```bash
# 解决方案：确保 CMake 能找到 pybind11
pip install pybind11
```

**问题：** `Python.h not found`
```bash
# 解决方案：安装 Python 开发包
sudo apt-get install python3-dev  # Ubuntu/Debian
```

### 运行时错误

**问题：** `ImportError: No module named '_radial_flow_cpp'`
```bash
# 解决方案：确保编译产物已复制到 python/ 目录
cp build/_radial_flow_cpp.*.so python/
```

**问题：** 计算不收敛
```python
# 解决方案：调整收敛参数
result = solver.solve(
    convergence_threshold=0.001,  # 更严格的阈值
    n_periods=5,                  # 更多检测周期
    max_r_factor=10.0             # 更大的积分范围
)
```

## 📚 参考文献

1. Pound & Wardell (2021), "Black hole perturbation theory and gravitational self-force"
2. 径向流方法原始论文（见项目中的 .tex 文件）

## 📝 许可证

本项目仅供学术研究使用。

## 👥 作者

- 开发者：[您的名字]
- 联系方式：[您的邮箱]

## 🙏 致谢

感谢 pybind11 项目提供的优秀 C++/Python 绑定工具。

---

**最后更新：** 2026-02-04
