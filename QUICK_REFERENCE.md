# 快速参考卡片

## 🚀 快速开始

```bash
# 1. 构建
./build.sh

# 2. 测试
python tests/test_simple_convergence.py

# 3. 使用
python scripts/solve_adaptive_high_freq.py
```

## 📊 基本使用

```python
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver

bh = SchwarzschildParams(r0=2.0)
wave = WaveParams(omega=0.5, l=2)
solver = RadialFlowSolver(bh, wave)

result = solver.solve(r_max=200.0)
print(f"|A/B| = {abs(result['reflection_coeff']):.6f}")
```

## ⚡ 高频计算（ω ≥ 3.0）

```python
# 方法 1：手动参数
result = solver.solve(
    r_max=200.0,
    epsilon=1e-3,      # ← 关键
    n_points=400000,   # ← 关键
    convergence_threshold=0.1,
    n_periods=10
)

# 方法 2：自适应（推荐）
from scripts.solve_adaptive_high_freq import solve_adaptive_high_freq
result = solve_adaptive_high_freq(bh, wave)
```

## 📋 参数速查表

| 频率 | epsilon | n_points | 说明 |
|-----|---------|----------|------|
| < 2.0 | 1e-6 | 10000 | 标准 |
| 2-3 | 1e-6 | 100000 | 增强 |
| 3-5 | 1e-3 | 400000 | 高频 |
| ≥ 5 | 1e-3 | 800000 | 极高频 |

## 🔍 调试命令

```bash
# 检查初值和系数
python tests/check_initial_values.py

# 测试高频
python tests/test_high_freq.py

# 步长收敛性
python tests/test_step_convergence.py

# 频率扫描
python scripts/frequency_sweep_adaptive.py
```

## 📚 文档速查

- **高频问题**：`docs/HIGH_FREQ_FINAL_SUMMARY.md`
- **收敛判据**：`docs/CONVERGENCE_FIX_SUMMARY.md`
- **完整报告**：`docs/FIX_COMPLETION_REPORT.md`
- **项目总览**：`README.md`

## 🎯 常见问题

### Q: 高频时出现 nan？
**A**: 增大 epsilon=1e-3，增加 n_points=400000

### Q: 收敛判据不工作？
**A**: 已修复！相位修正已应用

### Q: 如何选择参数？
**A**: 使用 `solve_adaptive_high_freq.py` 自动选择

## 📁 项目结构

```
radial_flow/
├── README.md          # 项目总览
├── docs/              # 📚 所有文档
├── scripts/           # 🔧 实用脚本
├── tests/             # 🧪 测试脚本
├── figures/           # 📊 图片
├── src/               # C++ 源码
├── include/           # C++ 头文件
└── python/            # Python 接口
```

## 💡 提示

- 从 `README.md` 开始
- 高频问题看 `docs/HIGH_FREQ_FINAL_SUMMARY.md`
- 使用 `scripts/solve_adaptive_high_freq.py` 自动处理
- 所有测试在 `tests/` 目录

---

**版本**: 1.0 | **状态**: ✅ 稳定 | **更新**: 2026-01-26
