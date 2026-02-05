# 文档索引

本目录包含径向流项目的所有文档。

## 📖 推荐阅读顺序

### 1. 新手入门

如果你是第一次使用这个项目：

1. **../README.md** - 项目总览和快速开始
2. **HIGH_FREQ_FINAL_SUMMARY.md** - 高频计算快速参考
3. **CONVERGENCE_FIX_SUMMARY.md** - 收敛判据说明

### 2. 深入理解

如果你想了解技术细节：

1. **FIX_COMPLETION_REPORT.md** - 完整修复报告
2. **HIGH_FREQ_CONVERGENCE_SOLUTION.md** - 高频收敛详细分析
3. **current_implementation.md** - 当前实现说明

### 3. 参考资料

如果你需要查阅特定信息：

1. **FREQUENCY_SWEEP_RESULTS.md** - 频率扫描结果
2. **J_and_sigma_definitions.md** - 数学定义
3. **PROJECT_DOCUMENTATION.md** - 项目文档

## 📚 文档分类

### 核心文档 ⭐

**HIGH_FREQ_FINAL_SUMMARY.md**
- 高频收敛问题的快速参考
- 推荐参数表
- 使用方法示例
- **适合**：需要快速解决高频计算问题

**CONVERGENCE_FIX_SUMMARY.md**
- 收敛判据修复总结
- 相位修正原理
- 修复前后对比
- **适合**：理解收敛检测如何工作

**FIX_COMPLETION_REPORT.md**
- 完整的修复报告
- 测试结果总结
- 物理验证
- **适合**：全面了解项目状态

### 技术文档 🔬

**HIGH_FREQ_CONVERGENCE_SOLUTION.md**
- 高频收敛问题的详细分析
- 刚性方程讨论
- 数值稳定性分析
- 测试结果对比
- **适合**：深入理解高频问题

**current_implementation.md**
- 当前代码实现说明
- 公式推导
- 代码对应关系
- **适合**：理解代码实现

**J_and_sigma_definitions.md**
- J 和 σ 的数学定义
- 从 .tex 文件提取
- **适合**：查阅数学公式

### 结果文档 📊

**FREQUENCY_SWEEP_RESULTS.md**
- 频率扫描详细结果
- 不同频率的行为分析
- 收敛特性讨论
- **适合**：了解物理结果

## 🎯 按问题查找

### 问题：高频时出现 nan

**查看**：HIGH_FREQ_FINAL_SUMMARY.md

**关键信息**：增大 epsilon（如 1e-3）+ 增加 n_points（如 400000-800000）

### 问题：收敛判据不工作

**查看**：CONVERGENCE_FIX_SUMMARY.md

**关键信息**：已修复！相位修正已应用

### 问题：如何选择参数

**查看**：HIGH_FREQ_FINAL_SUMMARY.md 或 ../README.md

**关键信息**：根据频率选择参数，或使用自适应求解器

---

**如有疑问，请先查看 HIGH_FREQ_FINAL_SUMMARY.md 和 CONVERGENCE_FIX_SUMMARY.md**
