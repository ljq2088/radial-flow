#!/bin/bash
# 一键编译脚本

echo "================================"
echo "Radial Flow 编译脚本"
echo "================================"

# 创建 build 目录
echo "创建 build 目录..."
mkdir -p build
cd build

# 运行 CMake
echo "运行 CMake 配置..."
cmake ..

# 编译
echo "编译 C++ 代码..."
make

# 复制编译产物到 python 目录
echo "复制编译产物到 python 目录..."
cp _radial_flow_cpp*.so ../python/

echo ""
echo "================================"
echo "编译完成！"
echo "================================"
echo ""
echo "现在可以运行演示脚本："
echo "  python3 demo.py"
echo "  python3 quick_test.py"
echo ""
