cmake_minimum_required(VERSION 3.29)
project(DPNN LANGUAGES CXX)

# -------------------- 编译选项 --------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -------------------- 编译优化选项 --------------------
# 设置构建类型为 Release（如果未指定）
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# 简化的优化选项
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3)
    add_compile_definitions(NDEBUG)
    add_compile_definitions(EIGEN_NO_DEBUG)
    message(STATUS "Release build with -O3 optimization")
endif()

# -------------------- 第三方库 --------------------
find_package(Boost REQUIRED)            # Boost库
find_package(CGAL REQUIRED)             # 提供 CGAL::CGAL
find_package(Eigen3 3.3 REQUIRED)       # 提供变量 EIGEN3_INCLUDE_DIR

# -------------------- 可执行文件 --------------------
add_executable(DPNN src/main.cpp src/dpnn.cpp)

# 设置包含目录
target_include_directories(DPNN
        PRIVATE
        ${EIGEN3_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# 链接库
target_link_libraries(DPNN
        PRIVATE
        CGAL::CGAL
        Boost::boost
)

# 启用链接时优化（如果支持）
set_property(TARGET DPNN PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

# -------------------- 额外信息 --------------------
message(STATUS "=== Build Configuration ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Using -O3 optimization for Release build")
message(STATUS "===============================")