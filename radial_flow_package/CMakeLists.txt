cmake_minimum_required(VERSION 3.15)
project(radial_flow LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找 pybind11
find_package(pybind11 REQUIRED)

# 包含目录
include_directories(${CMAKE_SOURCE_DIR}/include)

# C++ 核心库
add_library(radial_flow_core STATIC
    src/radial_flow.cpp
)

# 设置 PIC 标志
set_target_properties(radial_flow_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Python 扩展模块
pybind11_add_module(_radial_flow_cpp
    src/bindings.cpp
)

target_link_libraries(_radial_flow_cpp PRIVATE radial_flow_core)

# 安装
install(TARGETS _radial_flow_cpp
    LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/python
)