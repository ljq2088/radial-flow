#!/bin/bash
# 一键构建脚本

set -e  # 遇到错误立即退出

echo "=========================================="
echo "径向流计算项目 - 构建脚本"
echo "=========================================="

# 检查依赖
echo "检查依赖..."
python3 -c "import pybind11" 2>/dev/null || {
    echo "错误：未找到 pybind11，正在安装..."
    pip install pybind11
}

python3 -c "import numpy" 2>/dev/null || {
    echo "错误：未找到 numpy，正在安装..."
    pip install numpy
}

python3 -c "import scipy" 2>/dev/null || {
    echo "错误：未找到 scipy，正在安装..."
    pip install scipy
}

python3 -c "import matplotlib" 2>/dev/null || {
    echo "错误：未找到 matplotlib，正在安装..."
    pip install matplotlib
}

# 创建构建目录
echo "创建构建目录..."
mkdir -p build
cd build

# CMake 配置
echo "CMake 配置..."
cmake ..

# 编译
echo "编译 C++ 核心..."
make -j$(nproc)

# 安装
echo "安装 Python 模块..."
make install

cd ..

echo "=========================================="
echo "构建完成！"
echo "=========================================="
echo ""
echo "运行示例："
echo "  cd /home/ljq/code/radial_flow"
echo "  python -m python.radial_flow"
echo ""