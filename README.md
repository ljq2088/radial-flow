# 径向流计算项目

Schwarzschild 背景下标量场（s=0）的径向流方法计算反射系数。

## 项目结构

```
radial_flow/
├── include/
│   └── radial_flow.hpp          # C++ 头文件
├── src/
│   ├── radial_flow.cpp          # C++ 核心实现
│   └── bindings.cpp             # pybind11 Python 绑定
├── python/
│   ├── __init__.py              # Python 包初始化
│   └── radial_flow.py           # Python 接口层
├── tests/
│   └── test_scalar.py           # 单元测试
├── CMakeLists.txt               # CMake 构建配置
├── build.sh                     # 一键构建脚本
├── analysis.md                  # 推导验证和技术文档
└── README.md                    # 本文件
```

## 技术架构

- **C++ 核心**：高性能数值计算（RK4 积分器）
- **pybind11**：C++ 和 Python 之间的桥梁
- **Python 接口**：友好的 API 和可视化

## 快速开始

### 1. 安装依赖

```bash
pip install pybind11 numpy scipy matplotlib
```

### 2. 构建项目

```bash
cd /home/ljq/code/radial_flow
./build.sh
```

### 3. 运行示例

```bash
python -m python.radial_flow
```

## 使用示例

```python
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

# 设置参数
bh = SchwarzschildParams(r0=2.0)  # M = 1
wave = WaveParams(omega=0.5, l=2)

# 创建求解器
solver = RadialFlowSolver(bh, wave)

# 计算反射系数
reflection_coeff = solver.compute_reflection_coefficient()
print(f"A/B = {reflection_coeff}")
```

## 物理背景

计算 Schwarzschild 黑洞背景下标量场的反射系数 A/B：

1. **Riccati 方程**：$\partial_r\sigma = a + b\sigma + c\sigma^2$
2. **边界条件**：视界处 $\sigma_0 = 0$（纯入射波）
3. **渐近行为**：$\sigma \to -A/B$（反射系数）

详见 `analysis.md` 获取完整推导。

## 参考文献

- 径向流(check) (1).tex
- Pound & Wardell (2021) - Black hole perturbation theory

## 许可证

MIT License