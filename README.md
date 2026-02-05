# 径向流（Radial Flow）项目

黑洞扰动理论中的径向流方法实现，用于计算标量场在 Schwarzschild 黑洞中的反射系数。

## 项目结构

```
radial_flow/
├── README.md                 # 本文件
├── CMakeLists.txt           # CMake 构建配置
├── build.sh                 # 构建脚本
│
├── include/                 # C++ 头文件
│   └── radial_flow.hpp
│
├── src/                     # C++ 源文件
│   ├── radial_flow.cpp     # 核心实现
│   └── bindings.cpp        # Python 绑定
│
├── python/                  # Python 模块
│   ├── __init__.py
│   └── radial_flow.py      # Python 接口
│
├── build/                   # 构建目录（自动生成）
│
├── docs/                    # 📚 文档
│   ├── README.md           # 文档索引（推荐从这里开始）
│   ├── HIGH_FREQ_FINAL_SUMMARY.md          # 高频收敛快速参考 ⭐
│   ├── HIGH_FREQ_CONVERGENCE_SOLUTION.md   # 高频收敛详细文档
│   ├── CONVERGENCE_FIX_SUMMARY.md          # 收敛判据修复总结
│   ├── FIX_COMPLETION_REPORT.md            # 完整修复报告
│   ├── FREQUENCY_SWEEP_RESULTS.md          # 频率扫描结果
│   ├── PROJECT_DOCUMENTATION.md            # 项目文档
│   ├── current_implementation.md           # 当前实现说明
│   ├── J_and_sigma_definitions.md          # J 和 σ 定义
│   └── analysis.md                         # 分析文档
│
├── scripts/                 # 🔧 实用脚本
│   ├── solve_adaptive_high_freq.py         # 自适应高频求解器 ⭐
│   └── frequency_sweep_adaptive.py         # 频率扫描脚本
│
├── tests/                   # 🧪 测试脚本
│   ├── test_simple_convergence.py          # 简单收敛测试
│   ├── test_convergence.py                 # 多频率收敛测试
│   ├── test_high_freq.py                   # 高频测试
│   ├── test_step_convergence.py            # 步长收敛性测试
│   ├── test_different_epsilon.py           # epsilon 参数测试
│   ├── check_initial_values.py             # 检查初值和系数 ⭐
│   ├── check_high_freq_sigma.py            # 检查高频 σ 行为
│   └── ...                                 # 其他测试脚本
│
└── figures/                 # 📊 图片
    ├── frequency_sweep_adaptive.png        # 频率扫描结果
    ├── high_freq_omega3.0_check.png        # 高频检查图
    └── ...                                 # 其他图片
```

## 快速开始

### 1. 构建项目

```bash
./build.sh
```

或手动构建：
```bash
mkdir -p build && cd build
cmake ..
make -j$(nproc)
make install
cd ..
```

### 2. 基本使用

```python
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

# 创建黑洞和波参数
bh = SchwarzschildParams(r0=2.0)
wave = WaveParams(omega=0.5, l=2)
solver = RadialFlowSolver(bh, wave)

# 求解
result = solver.solve(r_max=200.0)

# 查看结果
print(f"|A/B| = {abs(result['reflection_coeff']):.6f}")
```

### 3. 高频计算（ω ≥ 3.0）

高频时需要特殊参数：

```python
# 方法 1：手动指定参数
result = solver.solve(
    r_max=200.0,
    epsilon=1e-3,      # 关键！
    n_points=400000,   # 关键！
    convergence_threshold=0.1,
    n_periods=10,
    max_r_factor=5.0
)

# 方法 2：使用自适应求解器（推荐）
from scripts.solve_adaptive_high_freq import solve_adaptive_high_freq
result = solve_adaptive_high_freq(bh, wave)
```

**为什么需要特殊参数？** 见 `docs/HIGH_FREQ_FINAL_SUMMARY.md`

## 核心功能

### 1. 反射系数计算

计算标量场在 Schwarzschild 黑洞中的反射系数 A/B：
- 低频（ω → 0）：|A/B| → 1（完全反射）
- 高频（ω → ∞）：|A/B| → 0（完全透射）

### 2. 自适应收敛

自动检测 σ 的收敛，无需手动指定 r_max：
- 基于相位修正后的震荡幅度
- 可配置阈值和检测周期数

### 3. 高频稳定性

通过增大 epsilon 和增加 n_points 处理刚性方程：
- ω < 2.0：标准参数
- 2.0 ≤ ω < 3.0：增强参数
- ω ≥ 3.0：高频参数

## 重要文档

### 📖 必读文档

1. **docs/HIGH_FREQ_FINAL_SUMMARY.md** ⭐
   - 高频收敛快速参考
   - 推荐参数表
   - 使用示例

2. **docs/CONVERGENCE_FIX_SUMMARY.md**
   - 收敛判据修复说明
   - 相位修正原理

3. **docs/FIX_COMPLETION_REPORT.md**
   - 完整修复报告
   - 测试结果总结

### 🔬 技术文档

- **docs/HIGH_FREQ_CONVERGENCE_SOLUTION.md** - 高频收敛详细分析
- **docs/current_implementation.md** - 当前实现说明
- **docs/J_and_sigma_definitions.md** - 数学定义

## 常用脚本

### 测试脚本

```bash
# 检查初值和 Riccati 系数（理解问题）
python tests/check_initial_values.py

# 简单收敛测试
python tests/test_simple_convergence.py

# 多频率测试
python tests/test_convergence.py

# 高频测试
python tests/test_high_freq.py

# 步长收敛性测试
python tests/test_step_convergence.py
```

### 实用脚本

```bash
# 自适应高频求解器（推荐）
python scripts/solve_adaptive_high_freq.py

# 频率扫描
python scripts/frequency_sweep_adaptive.py
```

## 调试指南

### 问题：高频时出现 nan

**症状**：ω ≥ 3.0 时，σ 变成 nan

**解决方案**：
1. 增大 epsilon（如 1e-3）
2. 增加 n_points（如 400000-800000）

**详细说明**：见 `docs/HIGH_FREQ_FINAL_SUMMARY.md`

### 问题：收敛判据不工作

**症状**：所有频率都显示"不收敛"

**解决方案**：已修复！相位修正已应用到收敛检测中。

**详细说明**：见 `docs/CONVERGENCE_FIX_SUMMARY.md`

## 参数说明

### solve() 参数

```python
result = solver.solve(
    r_max=200.0,              # 目标最大半径（单位：r0）
    epsilon=1e-6,             # 初始点离视界的距离
    n_points=10000,           # 积分点数
    convergence_threshold=0.01,  # 收敛阈值（1%）
    n_periods=2,              # 检测周期数
    max_r_factor=5.0          # 最大扩展因子
)
```

### 推荐参数

| 频率范围 | epsilon | n_points | threshold | n_periods |
|---------|---------|----------|-----------|-----------|
| ω < 2.0 | 1e-6 | 10000 | 0.01 | 2 |
| 2.0 ≤ ω < 3.0 | 1e-6 | 100000 | 0.01 | 5 |
| 3.0 ≤ ω < 5.0 | 1e-3 | 400000 | 0.1 | 10 |
| ω ≥ 5.0 | 1e-3 | 800000 | 0.1 | 10 |

## 物理结果

### 反射系数行为

- **低频**（ω = 0.1）：|A/B| = 1.000（完全反射）
- **中频**（ω = 0.5）：|A/B| = 0.767（部分反射）
- **高频**（ω = 3.0）：|A/B| = 0.027（主要透射）
- **极高频**（ω = 5.0）：|A/B| = 2.5×10⁻⁴（几乎完全透射）

符合物理预期：低频被势垒阻挡，高频穿透势垒。

## 技术细节

### 数值方法

- **积分方法**：RK4（4阶 Runge-Kutta）
- **边界条件**：视界处 σ₀ = 0（洛必达法则）
- **收敛检测**：基于相位修正后的震荡幅度

### 刚性方程处理

高频时 Riccati 方程变得刚性（stiff）：
- 系数 b ∝ ω，高频时 |b| ~ 10⁷
- RK4 需要非常小的步长才稳定
- 解决方案：增大 epsilon + 增加 n_points

详见：`docs/HIGH_FREQ_CONVERGENCE_SOLUTION.md`

## 贡献者

- 径向流方法实现
- 收敛判据修复
- 高频稳定性改进

## 许可证

（待添加）

## 参考文献

- 径向流(check) (1).tex - 理论公式来源
- Pound 和 Wardell 相关文献

---

**最后更新**：2026-01-26

**版本**：1.0

**状态**：✅ 收敛判据已修复，高频计算已稳定
