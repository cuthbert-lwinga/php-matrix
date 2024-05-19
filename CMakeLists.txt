cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Add Eigen directory
include_directories(${CMAKE_SOURCE_DIR}/eigen)

# Add your source files
add_executable(MyProject main.cpp matrix.cpp)