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
- ✅ **可移植**：使用相对路径，可在任何电脑上运行

## 📁 文件结构

```
radial_flow_package/
├── include/
│   └── radial_flow.hpp          # C++ 头文件
├── src/
│   ├── radial_flow.cpp          # C++ 实现
│   └── bindings.cpp             # Python 绑定
├── python/
│   └── radial_flow.py           # Python 包装层
├── build/                       # 编译输出目录（自动创建）
├── figures/                     # 图表输出目录
├── CMakeLists.txt               # CMake 配置文件
├── build.sh                     # 一键编译脚本
├── demo.py                      # 主演示脚本 ⭐
├── quick_test.py                # 快速测试脚本
├── example_usage.py             # 使用示例
├── .gitignore                   # Git 忽略文件
└── README.md                    # 本文件
```

## 🔧 依赖项

### 系统要求
- Linux / macOS / WSL
- C++17 或更高版本
- CMake 3.15+
- Python 3.8+

### Python 依赖
```bash
pip install numpy matplotlib
```

## 🚀 快速开始

### 步骤 1: 编译

```bash
./build.sh
```

### 步骤 2: 运行演示

```bash
python3 demo.py
```

### 步骤 3: 快速测试

```bash
python3 quick_test.py
```

## 📖 使用方法

### 基本用法

```python
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

# 设置参数
bh = SchwarzschildParams(r0=2.0)
wave = WaveParams(omega=0.2, l=0, use_subleading_term=True)

# 计算
solver = RadialFlowSolver(bh, wave)
result = solver.solve()

# 获取结果
print(f"|A/B| = {abs(result['reflection_coeff']):.6f}")
```

## 📊 示例输出

```
【简化版本】V = l(l+1)/r²
  |A/B| = 0.304564

【完整版本】V = l(l+1)/r² + f'/r
  |A/B| = 0.234955

【对比分析】
  相对差异: -22.86%
  物理意义: 次领头项使势垒增高，反射系数减小
```

## 📦 移植到其他电脑

1. 将整个 `radial_flow_package` 文件夹复制到目标电脑
2. 确保安装了依赖项
3. 运行 `./build.sh` 编译
4. 运行 `python3 demo.py` 测试

所有路径都使用相对路径，无需修改任何配置！

---

**最后更新：** 2026-02-05
